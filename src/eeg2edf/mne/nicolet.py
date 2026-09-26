"""Nicolet / Nervus .e as an mne.io.Raw."""
import numpy as np

from ..nicolet import convert, nicolet
from ._base import (EEG2EDFRaw, annotations, finish, guess_type, library_errors,
                    make_info, mne_name, orig_units, resolve_types, unique)


def _read_digital(ex, fname, start, stop):
    header, sfreq, n = ex["header"], ex["sfreq"], stop - start
    out = np.empty((len(ex["channels"]), n), dtype=np.int32)
    with open(fname, "rb") as fh:
        for j, (ch, rate, first) in enumerate(zip(ex["channels"], ex["rates"], ex["first"])):
            if rate == sfreq:
                lo = first + start
                v = np.frombuffer(nicolet.read_channel_range(fh, header, ch, lo, lo + n),
                                  dtype="<i2")
            else:
                # A slower channel is held: frame t shows the last sample at or before t.
                t = ex["offset_s"] + (start + np.arange(n)) / sfreq
                idx = np.minimum(np.floor(t * rate + 1e-9).astype(np.int64),
                                 header["channels"][ch]["n_samples"] - 1)
                lo, hi = int(idx[0]), int(idx[-1]) + 1
                v = np.frombuffer(nicolet.read_channel_range(fh, header, ch, lo, hi),
                                  dtype="<i2")[idx - lo]
            out[j] = convert.corrected(v, ex["invert"])
    return out


class RawNicolet(EEG2EDFRaw):
    """A Nicolet recording: the whole stored stream, or one segment.

    Use `read_raw_nicolet`.
    """

    def __init__(self, fname, segment=None, *, invert=True, include_trends=False,
                 ch_types=None, preload=False, verbose=None):
        fname = str(fname)
        with library_errors():
            header = nicolet.read_header(fname)
            channels = header["channels"]
            keep = convert.select_channels(channels, include_trends)
            pairs = convert.referential_pairs(channels, keep)
            segments = header["segments"]
            starts, at = [], 0.0
            for seg in segments:
                starts.append(at)
                at += seg["duration"]
            if segment is None:
                clip = (None, segments[0]["start"], 0.0, at)
            elif 0 <= segment < len(segments):
                clip = (segment, segments[segment]["start"], starts[segment],
                        segments[segment]["duration"])
            else:
                raise ValueError(f"segment {segment} out of range: {fname} has "
                                 f"{len(segments)}")
            meta = convert.build_sidecar(header, pairs, clip, invert, False)
        _index, _start, offset_s, duration = clip

        # One rate for the Raw: the fastest channel's. Slower ones (the 1 Hz
        # trends) are held sample-to-sample, which is what they mean anyway.
        sfreq = max(channels[a]["sfreq"] for _l, a, _b in pairs)
        n_times = int(round(duration * sfreq))
        units, scales, guessed, rates, first = [], [], [], [], []
        for label, a, _b in pairs:
            c = channels[a]
            unit = convert.EEG_UNIT if c["reference"] else convert.DERIVED_UNIT
            units.append(unit)
            scales.append(c["resolution"])
            guessed.append(guess_type(label, unit))
            rates.append(c["sfreq"])
            first.append(int(round(offset_s * c["sfreq"])))

        ch_names = unique([mne_name(label) for label, _a, _b in pairs])
        types = resolve_types(ch_names, guessed, ch_types)
        for c, n in zip(meta["channels"], ch_names):
            c["ch_name"] = n
        info = make_info(ch_names, sfreq, types, units, scales, meta)
        extras = {
            "header": header,
            "channels": [a for _l, a, _b in pairs],
            "rates": rates,
            "first": first,
            "offset_s": offset_s,
            "sfreq": sfreq,
            "invert": invert,
            "offsets": np.zeros(len(pairs)),
            "reader": _read_digital,
        }
        super().__init__(info, preload=preload, last_samps=(n_times - 1,),
                         filenames=[fname], raw_extras=[extras], orig_format="short",
                         orig_units=orig_units(ch_names, units, types), verbose=verbose)

        events = header["events"]
        # A label-less event of an unknown type carries nothing a reader could
        # use; like the EDF, the annotations leave it to the sidecar.
        def keep_event(i, _ev):
            return bool(events[i]["label"]) or not events[i]["type"].startswith("{")

        joins = starts[1:] if segment is None else []
        self.set_annotations(annotations(meta, ch_names, info["meas_date"],
                                         keep=keep_event, boundaries=joins))


def read_raw_nicolet(fname, segment=None, *, invert=True, include_trends=False,
                     montage=None, ch_types=None, positions=None, preload=False,
                     verbose=None):
    """Read a Nicolet / Nervus .e as an `mne.io.Raw`.

    Parameters
    ----------
    fname : path-like
        The .e file.
    segment : int | None
        One recorded segment, or None for the whole stored stream: the segments
        back to back, as the vendor exporter produced, with "BAD boundary" and
        "EDGE boundary" annotations at each join. Every segment's wall-clock
        start is in the sidecar's ``segments``.
    invert : bool
        The stored samples are inverted relative to uV; True (the default)
        corrects that, as nicolet2edf does.
    include_trends : bool
        Also read the derived trend channels (Rate, IBI, ...). They run at a
        lower rate and are held sample-to-sample onto the EEG rate, typed misc.
    montage : str | None
        The recording's display montage by name, returned instead of the
        referential channels. See `apply_montage`.
    ch_types, positions, preload
        As for `read_raw_nihon_kohden`.
    """
    raw = RawNicolet(fname, segment, invert=invert, include_trends=include_trends,
                     ch_types=ch_types, preload=preload, verbose=verbose)
    return finish(raw, montage=montage, positions=positions)
