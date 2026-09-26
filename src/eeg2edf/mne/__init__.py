"""MNE-Python readers: open a recording as an `mne.io.Raw`, no EDF in between.

    import eeg2edf
    raw = eeg2edf.read_raw("FILE.EEG")          # or .e, .vwr
    meta = eeg2edf.mne.get_sidecar(raw)          # the eeg2edf-sidecar/1 dict

Needs the `mne` extra: ``pip install "eeg2edf[mne]"``.

What MNE has a place for goes there: measurement date, sex and birth date,
device, filter settings, line frequency, channel types and calibration, and
every event as an annotation (with its sidecar fields as annotation extras).
Everything else the sidecar holds -- the vendor montages, segments,
per-channel filters and references, age -- is the sidecar itself, stored as
JSON in ``info["description"]`` so it survives a FIF round trip.
"""
import os

try:
    import mne  # noqa: F401
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        'eeg2edf.mne needs MNE-Python: pip install "eeg2edf[mne]"'
    ) from exc

from ._base import get_sidecar
from .micromed import RawMicromed, read_raw_micromed
from .montage import apply_montage, montage_names
from .nicolet import RawNicolet, read_raw_nicolet
from .nihon_kohden import RawNihonKohden, read_raw_nihon_kohden

READERS = {
    ".eeg": read_raw_nihon_kohden,
    ".e": read_raw_nicolet,
    ".vwr": read_raw_micromed,
}


def read_raw(fname, **kwargs):
    """Read any supported recording, choosing the reader by extension:
    .EEG `read_raw_nihon_kohden`, .e `read_raw_nicolet`, .vwr
    `read_raw_micromed`. Keyword arguments go to that reader."""
    ext = os.path.splitext(str(fname))[1].lower()
    if ext not in READERS:
        raise ValueError(f"don't know how to read {ext or 'a file with no extension'!r} "
                         f"files -- expected one of {', '.join(sorted(READERS))}")
    return READERS[ext](fname, **kwargs)


__all__ = [
    "RawMicromed", "RawNicolet", "RawNihonKohden", "apply_montage", "get_sidecar",
    "montage_names", "read_raw", "read_raw_micromed", "read_raw_nicolet",
    "read_raw_nihon_kohden",
]
