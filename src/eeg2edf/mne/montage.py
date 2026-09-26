"""Vendor display montages, applied to a Raw read by eeg2edf."""
import copy

import mne
import numpy as np

from .. import edfcommon
from ._base import get_sidecar, set_sidecar, unique


def montage_names(inst):
    """Names of the display montages the source recording carries."""
    return [m["name"] for m in get_sidecar(inst)["montages"]]


def _derivations(meta, montage, ch_names):
    """(label, active, reference, sign) per trace, as Raw channel names.

    Indices in the sidecar are the format's own (.21E codes, channel-table
    positions, LABCOD numbers) and so are each channel's `source_index`, which
    is what joins the two. An electrode the recording does not carry drops out
    of its trace the way the EDF converter drops it: Nihon Kohden's 0V
    reference contributes zero (`0V-EKG1` is -EKG1); Nicolet and Micromed
    traces need their active electrode and keep it alone when the reference
    is missing.
    """
    by_index = {}
    for c in meta["channels"]:
        name = c.get("ch_name")
        if not c["derived"] and c["source_index"] is not None and name in ch_names:
            by_index.setdefault(c["source_index"], name)
    allow_missing_active = meta["source"]["format"] == "nihon-kohden"
    out = []
    for t in montage["channels"]:
        a = by_index.get(t["active_index"])
        b = by_index.get(t["reference_index"])
        if a is None and (b is None or not allow_missing_active):
            continue
        out.append((t["label"], a, b))
    return out


def apply_montage(raw, name):
    """A new Raw holding a vendor montage's traces instead of the referential
    channels -- what the vendor viewer shows, and what `--montage` writes to EDF.

    `name` is one of `montage_names(raw)`. Stimulus and trigger channels are
    kept after the traces. The result is computed in float64 from the full
    data, so unlike the EDF output it is not quantised again. Annotations,
    measurement info and the sidecar (with `montage_applied` and `channels`
    updated) carry over.
    """
    meta = get_sidecar(raw)
    montage = next((m for m in meta["montages"] if m["name"] == name), None)
    if montage is None:
        raise ValueError(f"no montage {name!r} -- available: "
                         + ", ".join(repr(m["name"]) for m in meta["montages"]))
    traces = _derivations(meta, montage, raw.ch_names)
    if not traces:
        raise ValueError(f"montage {name!r} names no channel this Raw holds")

    labels = unique([label for label, _a, _b in traces])

    pos = {n: i for i, n in enumerate(raw.ch_names)}
    data = raw.get_data()
    rows = []
    for _label, a, b in traces:
        row = data[pos[a]].copy() if a is not None else np.zeros(data.shape[1])
        if b is not None:
            row -= data[pos[b]]
        rows.append(row)

    kept = [n for n, t in zip(raw.ch_names, raw.get_channel_types())
            if t in ("stim", "misc") and n not in labels]
    rows += [data[pos[n]] for n in kept]
    all_types = dict(zip(raw.ch_names, raw.get_channel_types()))
    types = [all_types[a if a is not None else b] for _l, a, b in traces]
    types += [all_types[n] for n in kept]

    info = mne.create_info(labels + kept, raw.info["sfreq"], types)
    src = raw.info
    with info._unlock():
        for key in ("meas_date", "subject_info", "device_info", "line_freq",
                    "highpass", "lowpass", "experimenter", "proj_name"):
            if src.get(key) is not None:
                info[key] = copy.deepcopy(src[key])
        for ch, n in zip(info["chs"][len(labels):], kept):
            old = src["chs"][pos[n]]
            ch.update(unit=old["unit"], cal=old["cal"], range=old["range"])
        for ch, (_l, a, b) in zip(info["chs"], traces):
            active = src["chs"][pos[a if a is not None else b]]
            ch["loc"] = active["loc"].copy()

    by_name = {c.get("ch_name"): c for c in meta["channels"]}
    channels = []
    for label, (_l, a, b) in zip(labels, traces):
        base = by_name.get(a if a is not None else b, {})
        channels.append(edfcommon._fill(edfcommon.CHANNEL, {
            "label": label, "edf_label": label, "ch_name": label,
            "source_index": base.get("source_index") if a is not None else None,
            "unit": base.get("unit"), "sfreq_hz": base.get("sfreq_hz"),
            "reference": b,
            "low_cut": base.get("low_cut"), "high_cut": base.get("high_cut"),
            "notch": base.get("notch"), "derived": b is not None,
        }))
    channels += [by_name[n] for n in kept if n in by_name]
    meta = dict(meta, channels=channels, montage_applied=montage["name"])
    set_sidecar(info, meta)

    out = mne.io.RawArray(np.asarray(rows), info, first_samp=raw.first_samp,
                          verbose="error")
    # An event pinned to a channel the montage consumed is no longer pinned;
    # its extras still name the channel.
    annot = raw.annotations
    names = set(out.ch_names)
    out.set_annotations(mne.Annotations(
        annot.onset, annot.duration, annot.description, orig_time=annot.orig_time,
        ch_names=[tuple(c for c in chs if c in names) for chs in annot.ch_names],
        extras=list(annot.extras)))
    return out

