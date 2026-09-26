"""Nihon Kohden EEG-1100 / EEG-1200A .EEG as an mne.io.Raw."""
import copy

import mne
import numpy as np

from ..nihon_kohden import convert, nk
from ._base import (EEG2EDFRaw, annotations, finish, get_sidecar, guess_type,
                    library_errors, make_info, mne_name, orig_units, resolve_types,
                    set_sidecar, unique)


def _scale(e21, dc_cal):
    """(unit, units per LSB, offset in units) for a channel word.

    The same calibration convert.signal_spec writes to the EDF header, as the
    exact per-LSB value rather than an 8-character physical range.
    """
    cal = dc_cal.get(convert.dc_key(e21) or "")
    if cal and cal.get("unit") and cal["coefficient"]:
        return cal["unit"], nk.DC_UV_PER_LSB / 1000 * cal["coefficient"], cal["offset"]
    if e21 in nk.DC_RANGE:
        return "mV", nk.DC_UV_PER_LSB / 1000, 0.0
    return "uV", nk.UV_PER_LSB, 0.0


def _read_digital(ex, fname, start, stop):
    frame, n = ex["frame"], stop - start
    with open(fname, "rb") as fh:
        fh.seek(ex["data_address"] + start * frame * 2)
        raw = np.fromfile(fh, dtype="<u2", count=n * frame)
    if raw.size != n * frame:
        raise EOFError(f"{fname}: short read at frame {start}")
    out = raw.reshape(n, frame)[:, ex["cols"]].T.astype(np.int32)
    # Channel words are offset-binary; the event word is a raw code,
    # reinterpreted as int16 unchanged (see nk.py).
    mark = ex["mark"]
    out[~mark] -= 32768
    out[mark] = np.where(out[mark] > 32767, out[mark] - 65536, out[mark])
    return out


class RawNihonKohden(EEG2EDFRaw):
    """One datablock (clip) of a Nihon Kohden recording.

    Use `read_raw_nihon_kohden`, which also joins every block with
    ``block="all"``.
    """

    def __init__(self, fname, block=0, *, all_channels=False, log=None, no_log=False,
                 annotations_csv=None, marker_channel=True, ch_types=None,
                 preload=False, verbose=None, _rec=None):
        fname = str(fname)
        with library_errors():
            rec = _rec or convert.load(fname, all_channels=all_channels, log=log,
                                       no_log=no_log, annotations=annotations_csv)
            blocks = rec["blocks"]
            if not 0 <= block < len(blocks):
                raise ValueError(f"block {block} out of range: {fname} has {len(blocks)}")
            blk, keep, dc_cal = blocks[block], rec["keep"], rec["dc_cal"]
            names, pairs = convert.signal_pairs(blk, keep)
            signals = convert.signal_specs(blk, keep, names, dc_cal, None, marker_channel)
            events = convert.block_events(rec, block, 1)
        meta = convert.build_sidecar(fname, blk, names, pairs, signals, rec["e21"],
                                     rec["patient"], rec["montages"], None, events)
        meta["clip"]["index"] = block
        meta["clip"]["duration_s"] = blk["n_samples"] / blk["sfreq"]

        units, scales, offsets, cols, mark, guessed = [], [], [], [], [], []
        for label, (a, _b) in zip(names, pairs):
            unit, scale, offset = _scale(blk["e21_index"][a], dc_cal)
            units.append(unit)
            scales.append(scale)
            offsets.append(offset / scale)
            cols.append(a)
            mark.append(False)
            guessed.append(guess_type(label, unit))
        labels = list(names)
        if marker_channel:
            labels.append(convert.EVENT_LABEL)
            units.append("")
            scales.append(1.0)
            offsets.append(0.0)
            cols.append(blk["n_channels"])
            mark.append(True)
            guessed.append("stim")

        ch_names = unique([mne_name(n) for n in labels])
        types = resolve_types(ch_names, guessed, ch_types)
        for c, n in zip(meta["channels"], ch_names):
            c["ch_name"] = n
        info = make_info(ch_names, blk["sfreq"], types, units, scales, meta)
        extras = {
            "data_address": blk["data_address"],
            "frame": blk["n_channels"] + 1,
            "cols": np.asarray(cols, dtype=np.intp),
            "mark": np.asarray(mark, dtype=bool),
            "offsets": np.asarray(offsets, dtype=np.float64),
            "reader": _read_digital,
        }
        super().__init__(info, preload=preload, last_samps=(blk["n_samples"] - 1,),
                         filenames=[fname], raw_extras=[extras], orig_format="short",
                         orig_units=orig_units(ch_names, units, types), verbose=verbose)
        self.set_annotations(annotations(meta, ch_names, info["meas_date"]))


def _join_sidecars(raws):
    """One sidecar for back-to-back blocks: each block becomes a segment and
    event onsets move onto the joined timeline."""
    metas = [get_sidecar(r) for r in raws]
    meta = copy.deepcopy(metas[0])
    sfreq = raws[0].info["sfreq"]
    segments, events, at = [], [], 0.0
    for i, (m, r) in enumerate(zip(metas, raws)):
        duration = r.n_times / sfreq
        segments.append({"index": m["clip"]["index"], "start": m["clip"]["start"],
                         "offset_s": at, "duration_s": duration})
        for ev in m["events"]:
            ev = dict(ev)
            if ev["onset_s"] is not None:
                ev["onset_s"] += at
            events.append(ev)
        at += duration
    meta["clip"].update(index=None, offset_s=0.0, duration_s=at)
    meta["segments"] = segments
    meta["events"] = events
    return meta


def read_raw_nihon_kohden(fname, block=0, *, all_channels=False, log=None, no_log=False,
                          annotations=None, marker_channel=True, montage=None,
                          ch_types=None, positions=None, preload=False, verbose=None):
    """Read a Nihon Kohden .EEG (EEG-1100 or EEG-1200A) as an `mne.io.Raw`.

    Parameters
    ----------
    fname : path-like
        The .EEG file. Its siblings sharing the stem are used when present:
        .21E (names, reference, device), .PNT (sex, birth date), .LOG and .sld
        (events), .11D (DC calibration), .PTN folder (montages).
    block : int | "all"
        Which datablock (clip) to read. "all" joins every block back to back,
        with "BAD boundary"/"EDGE boundary" annotations at the joins; each
        block's true start stays in the sidecar's ``segments``.
    all_channels : bool
        Keep inputs the recording marks as unconnected, as ``--all-channels``.
    log : path-like | None
        A .LOG to use instead of the one next to the .EEG.
    no_log : bool
        Ignore the .LOG.
    annotations : path-like | None
        Extra events, a CSV of ``<ISO datetime|seconds>,<text>`` (see nk2edf).
    marker_channel : bool
        Keep the per-frame event word as a stim channel "Events/Markers".
    montage : str | None
        A .PTN montage to return instead of referential channels, or "auto" for
        the one the .LOG names at the block's start. See `apply_montage`.
    ch_types : str | dict | None
        Override the guessed channel types: one type (e.g. "seeg") for every
        channel guessed as EEG, or a dict of channel name to type.
    positions : str | DigMontage | None
        Sensor positions, e.g. "standard_1020". Off by default: SEEG contact
        names such as "A1" collide with 10-20 electrode names.
    preload : bool
        Load the data now. Otherwise it is read from disk on demand.
    """
    with library_errors():
        rec = convert.load(str(fname), all_channels=all_channels, log=log,
                           no_log=no_log, annotations=annotations)
        kw = dict(all_channels=all_channels, marker_channel=marker_channel,
                  ch_types=ch_types, preload=preload, verbose=verbose, _rec=rec)
        if block == "all":
            if len({b["sfreq"] for b in rec["blocks"]}) > 1:
                raise ValueError("blocks differ in sampling rate; read them one at a time")
            if len({tuple(b["e21_index"]) for b in rec["blocks"]}) > 1:
                raise ValueError("blocks differ in channels; read them one at a time")
            if any(e.get("seconds") is not None for e in rec["log_events"]):
                raise ValueError("annotations given in seconds are ambiguous across "
                                 "blocks; use ISO datetimes or read one block")
            raws = [RawNihonKohden(fname, i, **kw) for i in range(len(rec["blocks"]))]
            meta = _join_sidecars(raws)
            raw = mne.concatenate_raws(raws, verbose="error")
            set_sidecar(raw.info, meta)
            first = 0
        else:
            raw = RawNihonKohden(fname, block, **kw)
            first = block
        if montage is not None and montage.lower() == "auto":
            montage = convert.resolve_montage(rec, first, "auto")["name"]
    return finish(raw, montage=montage, positions=positions)
