"""Readers and EDF+ converters for proprietary clinical EEG formats.

The converters are the `nk2edf`, `nicolet2edf` and `vwr2edf` commands (or
`eeg2edf`, which picks one by file extension). With the `mne` extra
installed, `read_raw` opens any of the formats as an `mne.io.Raw`.
"""

__version__ = "0.1.0"


def read_raw(fname, **kwargs):
    """Open a recording as an `mne.io.Raw`; see `eeg2edf.mne.read_raw`."""
    from .mne import read_raw as _read_raw

    return _read_raw(fname, **kwargs)
