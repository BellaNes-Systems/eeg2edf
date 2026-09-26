"""Micromed VWR / Brain-Quick .vwr as an mne.io.Raw."""
import numpy as np

from ..micromed import convert, vwr
from ._base import (EEG2EDFRaw, annotations, finish, guess_type, library_errors,
                    make_info, mne_name, orig_format, orig_units, resolve_types,
                    unique)


def _read_digital(ex, _fname, start, stop):
    return vwr.read_frames(ex["header"], start, stop).T.astype(np.int64)


class RawMicromed(EEG2EDFRaw):
    """A Micromed VWR recording. Use `read_raw_micromed`."""

    def __init__(self, fname, *, ch_types=None, preload=False, verbose=None):
        fname = str(fname)
        with library_errors():
            header = vwr.read_header(fname)
            pairs = convert.referential_pairs(header)
            edf_labels = convert._unique_labels([label for label, _a, _b in pairs])
            low, high = convert._signal_ranges(header, pairs)
            meta = convert.build_sidecar(header, pairs, edf_labels, low, high, None)
        meta["clip"]["index"] = None
        meta["clip"]["offset_s"] = 0.0

        # Channels are named by electrode; each one's recording ground stays
        # in the sidecar as its reference. Labels that only the ground tells
        # apart keep the ground in the name.
        labels = [c.label for c in header.channels]
        if len(set(labels)) != len(labels):
            labels = [c.label if labels.count(c.label) == 1 else c.name
                      for c in header.channels]
        ch_names = unique([mne_name(n or f"CH{i + 1}") for i, n in enumerate(labels)])
        units = ["uV"] * len(ch_names)
        types = resolve_types(ch_names, [guess_type(n) for n in ch_names], ch_types)
        for c, n in zip(meta["channels"], ch_names):
            c["ch_name"] = n
        info = make_info(ch_names, header.frequency, types, units,
                         [c.factor for c in header.channels], meta)
        extras = {
            "header": header,
            "offsets": np.asarray([-c.lground for c in header.channels], dtype=np.float64),
            "reader": _read_digital,
        }
        super().__init__(info, preload=preload, last_samps=(header.n_samples - 1,),
                         filenames=[fname], raw_extras=[extras],
                         orig_format=orig_format(header.int_size),
                         orig_units=orig_units(ch_names, units, types), verbose=verbose)

        end = header.n_samples / header.frequency
        # TRONCA parts are where acquisition was cut and resumed.
        joins = [p.frame / header.frequency for p in header.parts
                 if 0 < p.frame < header.n_samples]
        self.set_annotations(annotations(
            meta, ch_names, info["meas_date"],
            keep=lambda _i, ev: ev["onset_s"] < end, boundaries=joins))


def read_raw_micromed(fname, *, montage=None, ch_types=None, positions=None,
                      preload=False, verbose=None):
    """Read a Micromed VWR / Brain-Quick .vwr as an `mne.io.Raw`.

    Parameters
    ----------
    fname : path-like
        The .vwr file.
    montage : str | None
        A stored montage by name, returned instead of the referential
        channels, or "auto" for the one HISTORY says was recorded. See
        `apply_montage`.
    ch_types, positions, preload
        As for `read_raw_nihon_kohden`.
    """
    raw = RawMicromed(fname, ch_types=ch_types, preload=preload, verbose=verbose)
    if montage is not None and montage.lower() == "auto":
        montage = raw._raw_extras[0]["header"].recorded_montage
        if montage is None:
            raise ValueError("no montage recorded in HISTORY; pass montage=NAME instead")
    return finish(raw, montage=montage, positions=positions)
