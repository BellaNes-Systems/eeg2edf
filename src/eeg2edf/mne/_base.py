"""What the three MNE readers share: the Raw base class, Info population from
a sidecar, channel typing, and the sidecar kept in info["description"]."""
import contextlib
import datetime as dt
import functools
import json
import re

import mne
import numpy as np
from mne._fiff.constants import FIFF
from mne._fiff.utils import _mult_cal_one

from .. import edfcommon

# Rows decoded per step when a read is large: every step builds a float64
# temporary of the whole step, so bounding it keeps preload memory at the size
# of the output rather than several times it.
CHUNK_FRAMES = 1 << 16

_UNIT_SCALE = {"uV": 1e-6, "µV": 1e-6, "mV": 1e-3, "V": 1.0}


@contextlib.contextmanager
def library_errors():
    """The converters stop with SystemExit, which is right for a command and
    wrong for a library call; readers raise ValueError instead."""
    try:
        yield
    except SystemExit as exc:
        raise ValueError(str(exc.code)) from None


class EEG2EDFRaw(mne.io.BaseRaw):
    """Base for the readers. A subclass puts a function in its raw extras,
    ``reader(raw_extras, filename, start, stop) -> ndarray``, returning the
    stored integers for frames [start, stop) of every channel in channel
    order, plus the per-channel ``offsets`` (in LSB) to add before the
    calibration is applied.

    It has to be a plain function: MNE reads through a proxy that exposes
    only ``_raw_extras`` and ``filenames``, not the instance.
    """

    def _read_segment_file(self, data, idx, fi, start, stop, cals, mult):
        extras, fname = self._raw_extras[fi], self._filenames[fi]
        read = extras["reader"]
        at = 0
        for lo in range(start, stop, CHUNK_FRAMES):
            hi = min(lo + CHUNK_FRAMES, stop)
            block = read(extras, fname, lo, hi).astype(np.float64)
            block += extras["offsets"][:, None]
            _mult_cal_one(data[:, at:at + hi - lo], block, idx, cals, mult)
            at += hi - lo


# --- the sidecar in info["description"] -------------------------------------

def get_sidecar(inst):
    """The `eeg2edf-sidecar/1` metadata of a Raw (or Info) read by eeg2edf.

    It lives in ``info["description"]`` as JSON, so it survives ``copy()``,
    cropping, picking and a FIF save/load round trip. It describes the channels
    as they were read: picking or renaming later does not rewrite it.
    """
    info = inst.info if hasattr(inst, "info") else inst
    text = info.get("description") or ""
    try:
        meta = json.loads(text)
    except ValueError:
        meta = None
    if not isinstance(meta, dict) or meta.get("schema") != edfcommon.SCHEMA:
        raise ValueError(
            "info['description'] does not hold an eeg2edf sidecar -- was this "
            "Raw read with eeg2edf, and has its description been overwritten?"
        )
    return meta


def set_sidecar(info, meta):
    info["description"] = json.dumps(meta, separators=(",", ":"), default=str)


# --- channels ----------------------------------------------------------------

def standard_montage_name(system="1005"):
    """MNE 1.13 renamed standard_10xx to colin27_10xx; use what this MNE has."""
    builtin = mne.channels.get_builtin_montages()
    new = f"colin27_{system}"
    return new if new in builtin else f"standard_{system}"


@functools.lru_cache(maxsize=1)
def _standard_names():
    names = mne.channels.make_standard_montage(standard_montage_name()).ch_names
    return {n.lower(): n for n in names}


def mne_name(label):
    """Nihon Kohden's built-in names carry an "EEG " prefix ("EEG FP1"). Drop
    it, in the standard spelling, where what is left is a 10-05 electrode --
    that is what lets set_montage find it. Anything else is left alone."""
    if label.upper().startswith("EEG "):
        std = _standard_names().get(label[4:].strip().lower())
        if std:
            return std
    return label


_ECG = re.compile(r"ECG|EKG")


def guess_type(label, unit="uV"):
    """MNE channel type from a label and physical unit."""
    u = label.upper()
    if _ECG.search(u):
        return "ecg"
    if "EOG" in u:
        return "eog"
    if "EMG" in u:
        return "emg"
    if unit not in ("uV", "µV"):
        return "misc"
    return "eeg"


def unique(names):
    out, seen = [], {}
    for n in names:
        if n in seen:
            seen[n] += 1
            n = f"{n}-{seen[n]}"
        else:
            seen[n] = 0
        out.append(n)
    return out


def resolve_types(names, guessed, ch_types):
    """Apply the user's `ch_types`: a dict of name -> type, or one type
    (e.g. "seeg") for every channel that was guessed as EEG."""
    types = list(guessed)
    if ch_types is None:
        return types
    if isinstance(ch_types, str):
        return [ch_types if t == "eeg" else t for t in types]
    unknown = set(ch_types) - set(names)
    if unknown:
        raise ValueError(f"ch_types names channels that are not in the recording: "
                         f"{sorted(unknown)}")
    return [ch_types.get(n, t) for n, t in zip(names, types)]


def make_info(names, sfreq, types, units, scales, meta):
    """Info for one recording. `units` are the physical units the scales are
    in ("uV", "mV", ... or "" for a dimensionless number); `scales` are
    physical units per stored LSB."""
    info = mne.create_info(names, float(sfreq), types)
    for ch, unit, scale, kind in zip(info["chs"], units, scales, types):
        to_si = _UNIT_SCALE.get(unit)
        ch["range"] = 1.0
        if kind == "stim":
            ch["cal"] = 1.0
            ch["unit"] = FIFF.FIFF_UNIT_NONE
        elif to_si is None:
            ch["cal"] = float(scale)
            ch["unit"] = FIFF.FIFF_UNIT_NONE
        else:
            ch["cal"] = float(scale) * to_si
            ch["unit"] = FIFF.FIFF_UNIT_V

    with info._unlock():
        start = meta["clip"]["start"]
        if start:
            # Vendor clocks are local wall time with no zone; MNE's EDF reader
            # makes the same choice of calling that UTC.
            info["meas_date"] = dt.datetime.fromisoformat(start).replace(
                tzinfo=dt.timezone.utc)
        subject = {}
        sex = (meta["patient"].get("sex") or "").upper()[:1]
        if sex in ("M", "F"):
            subject["sex"] = {"M": FIFF.FIFFV_SUBJ_SEX_MALE,
                              "F": FIFF.FIFFV_SUBJ_SEX_FEMALE}[sex]
        if meta["patient"].get("dob"):
            subject["birthday"] = dt.date.fromisoformat(meta["patient"]["dob"])
        if subject:  # never a name or record number
            info["subject_info"] = subject
        info["device_info"] = {"type": meta["source"]["format"]}
        if meta.get("device"):
            info["device_info"]["model"] = str(meta["device"])

        eeg = [c for c, t in zip(meta["channels"], types) if t in ("eeg", "seeg", "ecog")]
        lows = [c["low_cut"] for c in eeg if c.get("low_cut")]
        highs = [c["high_cut"] for c in eeg if c.get("high_cut")]
        if lows and len(lows) == len(eeg):
            info["highpass"] = float(max(lows))
        if highs and len(highs) == len(eeg):
            info["lowpass"] = float(min(highs))
        notches = {c.get("notch") for c in eeg}
        if len(notches) == 1 and notches <= {50, 60}:
            info["line_freq"] = float(notches.pop())
    set_sidecar(info, meta)
    return info


def orig_units(names, units, types):
    return {n: ("n/a" if t == "stim" or not u else u) for n, u, t in zip(names, units, types)}


def orig_format(nbytes):
    return "short" if nbytes <= 2 else "int"


# --- events ------------------------------------------------------------------

_ANNOT_EXTRA_TYPES = (int, float, str, type(None))
_ANNOT_RESERVED = {"onset", "duration", "description", "ch_names", "orig_time"}


def annotations(meta, ch_names, orig_time, keep=None, boundaries=()):
    """Annotations for every sidecar event that lands in the clip, plus MNE's
    "BAD boundary"/"EDGE boundary" pair at each join in the stored stream.

    `keep(i, event)` can veto an event. Sidecar fields other than onset,
    duration and label ride along as annotation extras.
    """
    onset, duration, description, chans, extras = [], [], [], [], []
    for i, ev in enumerate(meta["events"]):
        if ev["onset_s"] is None or (keep and not keep(i, ev)):
            continue
        onset.append(ev["onset_s"])
        duration.append(ev["duration_s"] or 0.0)
        description.append(ev["label"] or ev["type"] or "")
        ch = ev.get("channel")
        chans.append((ch,) if ch and ch in ch_names else ())
        extras.append({
            k: v for k, v in ev.items()
            if k not in ("onset_s", "duration_s", "label") and k not in _ANNOT_RESERVED
            and isinstance(v, _ANNOT_EXTRA_TYPES) and not isinstance(v, bool)
        })
    for at in boundaries:
        for text in ("BAD boundary", "EDGE boundary"):
            onset.append(at)
            duration.append(0.0)
            description.append(text)
            chans.append(())
            extras.append({})
    return mne.Annotations(onset, duration, description, orig_time=orig_time,
                           ch_names=chans, extras=extras)


def finish(raw, *, montage=None, positions=None):
    """What every reader does last: sensor positions, then a vendor montage."""
    if positions is not None:
        if isinstance(positions, str):
            positions = mne.channels.make_standard_montage(positions)
        raw.set_montage(positions, match_case=False, on_missing="ignore")
    if montage is not None:
        from .montage import apply_montage

        raw = apply_montage(raw, montage)
    return raw
