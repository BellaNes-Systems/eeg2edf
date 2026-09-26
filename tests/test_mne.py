"""The MNE readers against the EDF converters, on synthetic recordings.

The converters are verified against vendor software (see docs/); these tests
hold the readers to the converters, so both paths give the same signal."""
import datetime as dt
import json
import subprocess
import sys

import numpy as np
import pytest

mne = pytest.importorskip("mne")

import synth  # noqa: E402
from test_sidecar_schema import TOP  # noqa: E402

import eeg2edf  # noqa: E402
from eeg2edf import edfcommon  # noqa: E402
from eeg2edf.micromed import convert as vwr2edf  # noqa: E402
from eeg2edf.mne import apply_montage, get_sidecar, montage_names  # noqa: E402
from eeg2edf.nicolet import convert as nicolet2edf  # noqa: E402
from eeg2edf.nicolet import nicolet  # noqa: E402
from eeg2edf.nihon_kohden import convert as nk2edf  # noqa: E402


@pytest.fixture(scope="module")
def files(tmp_path_factory):
    d = tmp_path_factory.mktemp("rec")
    nk = synth.write_nk(d)
    nic = synth.write_nicolet(d)
    synth.write_vwr(d / "SYN.vwr")
    return {"nk": nk["eeg"], "nk_frames": nk["blocks"], "nic": nic["path"],
            "nic_streams": nic["streams"], "vwr": d / "SYN.vwr", "dir": d}


def edf(path):
    return mne.io.read_raw_edf(path, preload=True, verbose="error")


def converted(tool, *args):
    tool.main([str(a) for a in args])


def assert_within(actual, expected, tol):
    """|actual - expected| <= tol, with tol per channel (a column)."""
    err = np.abs(np.asarray(actual) - np.asarray(expected))
    worst = (err - tol).max()
    assert worst <= 1e-12, f"exceeds tolerance by {worst:g}; max error {err.max():g}"


SI = {"uV": 1e-6, "mV": 1e-3}


def lsb(meta):
    """Per-channel volts per stored LSB, from the sidecar. The resolution is in
    the channel's own unit (mV for a Nihon Kohden DC input)."""
    return np.array([c["resolution_uv_per_lsb"] * SI.get(c["unit"], 1.0)
                     for c in meta["channels"]])


# --- same signal as the EDF converters -------------------------------------

def test_nihon_kohden_matches_nk2edf(files, tmp_path):
    converted(nk2edf, files["nk"], tmp_path, "--blocks", "1")
    ref = edf(next(tmp_path.glob("*.edf")))
    raw = eeg2edf.read_raw(files["nk"], block=1)
    assert raw.ch_names == ref.ch_names
    # Events/Markers is a stim channel here: raw codes, where the EDF spreads
    # them over physical -1..1 (tested against the stored word below).
    tol = lsb(get_sidecar(raw))[:-1, None] / 2
    assert_within(raw.get_data()[:-1], ref.get_data()[:-1], tol)


def test_nihon_kohden_values_are_the_stored_words(files):
    raw = eeg2edf.read_raw(files["nk"], block=0)
    frames = files["nk_frames"][0].astype(np.int64)
    np.testing.assert_allclose(raw.get_data(picks="Fp1")[0] * 1e6,
                               (frames[:, 0] - 32768) * 6400 / 65536)
    marks = raw.get_data(picks="Events/Markers")[0]
    assert marks[10] == 1 and marks[20] == -1 and np.count_nonzero(marks) == 2
    assert raw.get_channel_types() == ["eeg", "eeg", "eeg", "ecg", "misc", "stim"]


def test_nicolet_matches_nicolet2edf_concat(files, tmp_path):
    converted(nicolet2edf, files["nic"], tmp_path, "--concat")
    ref = edf(tmp_path / "SAMPLE.edf")
    raw = eeg2edf.read_raw(files["nic"])
    assert raw.ch_names == ref.ch_names
    np.testing.assert_allclose(raw.get_data(), ref.get_data(), rtol=0, atol=1e-12)
    # inverted, with -32768 clamped rather than wrapped
    stored = files["nic_streams"][0].astype(np.int64)
    np.testing.assert_allclose(raw.get_data(picks="Fp1")[0] * 1e6,
                               -np.clip(stored, -32767, None) * 0.25)


def test_nicolet_segment_is_a_slice_of_the_stream(files):
    whole = eeg2edf.read_raw(files["nic"]).get_data()
    second = eeg2edf.read_raw(files["nic"], segment=1)
    np.testing.assert_array_equal(second.get_data(), whole[:, 512:768])
    assert second.info["meas_date"] == dt.datetime(2022, 9, 21, 3, 50, 41,
                                                   tzinfo=dt.timezone.utc)
    assert not any("boundary" in d for d in second.annotations.description)


def test_nicolet_trends_are_held_onto_the_eeg_rate(files):
    raw = eeg2edf.read_raw(files["nic"], include_trends=True)
    assert raw.ch_names[-1] == "Rate"
    assert raw.get_channel_types()[-1] == "misc"
    rate = raw.get_data(picks="Rate")[0]
    stored = -np.clip(files["nic_streams"][3].astype(np.int64), -32767, None)
    np.testing.assert_array_equal(rate, np.repeat(stored, 256))


def test_micromed_matches_vwr2edf(files, tmp_path):
    converted(vwr2edf, files["vwr"], tmp_path)
    ref = edf(tmp_path / "SYN.edf")
    raw = eeg2edf.read_raw(files["vwr"])
    assert raw.ch_names == ["B", "A"]  # by electrode; the ground is the reference
    assert [c["reference"] for c in get_sidecar(raw)["channels"]] == ["REF", "REF"]
    np.testing.assert_allclose(raw.get_data(), ref.get_data(), atol=1e-9)
    np.testing.assert_allclose(raw.get_data() * 1e6, [[1, 2, 3, 4], [1, 2, 3, 4]])


# --- lazy reading ----------------------------------------------------------

@pytest.mark.parametrize("key", ["nk", "nic", "vwr"])
def test_lazy_equals_preload(files, key, monkeypatch):
    from eeg2edf.mne import _base

    full = eeg2edf.read_raw(files[key], preload=True).get_data()
    monkeypatch.setattr(_base, "CHUNK_FRAMES", 3)  # force many small reads
    lazy = eeg2edf.read_raw(files[key])
    np.testing.assert_array_equal(lazy.get_data(), full)
    tmax = (lazy.n_times - 2) / lazy.info["sfreq"]
    cropped = lazy.copy().crop(1 / lazy.info["sfreq"], tmax)
    np.testing.assert_array_equal(cropped.get_data(), full[:, 1:-1])
    picked = eeg2edf.read_raw(files[key]).pick([lazy.ch_names[-1]])
    np.testing.assert_array_equal(picked.get_data()[0], full[-1])


def test_read_channel_range_follows_the_index(files):
    header = nicolet.read_header(files["nic"])
    with open(files["nic"], "rb") as fh:
        whole = nicolet.read_channel(fh, header, 0)
        for lo, hi in ((0, 768), (0, 1), (300, 500), (383, 385), (767, 768), (5, 5)):
            assert nicolet.read_channel_range(fh, header, 0, lo, hi) == whole[lo * 2:hi * 2]


# --- metadata --------------------------------------------------------------

def test_nihon_kohden_info(files):
    raw = eeg2edf.read_raw(files["nk"])
    info = raw.info
    assert info["meas_date"] == dt.datetime(2024, 2, 6, 19, 55, 47, tzinfo=dt.timezone.utc)
    assert info["subject_info"]["sex"] == 1
    assert info["subject_info"]["birthday"] == dt.date(2010, 5, 1)
    assert "his_id" not in info["subject_info"] and "last_name" not in info["subject_info"]
    assert info["device_info"] == {"type": "nihon-kohden", "model": "EEG-1200A"}
    meta = get_sidecar(raw)
    assert meta["patient"]["age_at_recording"] == "13"
    assert meta["reference"] == "C3C4"
    assert montage_names(raw) == ["EMU1"]


def test_nicolet_filters_reach_info(files):
    info = eeg2edf.read_raw(files["nic"]).info
    assert (info["highpass"], info["lowpass"], info["line_freq"]) == (0.5, 70.0, 60.0)


@pytest.mark.parametrize("key", ["nk", "nic", "vwr"])
def test_sidecar_is_the_converter_schema(files, key):
    raw = eeg2edf.read_raw(files[key])
    meta = get_sidecar(raw)
    assert list(meta)[:len(TOP)] == TOP
    for section, template in (("channels", edfcommon.CHANNEL), ("events", edfcommon.EVENT),
                              ("segments", edfcommon.SEGMENT)):
        for row in meta[section]:
            assert list(row)[:len(template)] == list(template)
    assert [c["ch_name"] for c in meta["channels"]] == raw.ch_names
    dt.datetime.fromisoformat(meta["clip"]["start"])


@pytest.mark.parametrize("key", ["nk", "nic", "vwr"])
def test_annotations_are_the_sidecar_events(files, key):
    raw = eeg2edf.read_raw(files[key])
    meta = get_sidecar(raw)
    annot = raw.annotations
    real = [i for i, d in enumerate(annot.description) if "boundary" not in d]
    placed = [e for e in meta["events"]
              if e["onset_s"] is not None and e["onset_s"] < raw.times[-1] + 1]
    assert sorted((annot.onset[i], annot.description[i]) for i in real) == \
        sorted((e["onset_s"], e["label"]) for e in placed)
    for i in real:
        assert annot.extras[i]["source"]


def test_nicolet_event_channel_and_gap(files):
    raw = eeg2edf.read_raw(files["nic"])
    annot = raw.annotations
    i = list(annot.description).index("Seizure - big one")
    assert annot.ch_names[i] == ("Fp1",)
    assert annot.extras[i]["guid"] == synth.NIC_SEIZURE
    assert "in gap" not in annot.description  # between segments: sidecar only
    gap = next(e for e in get_sidecar(raw)["events"] if e["label"] == "in gap")
    assert gap["onset_s"] is None and gap["when"] == "2022-09-21T03:50:21"
    assert list(annot.onset[[d.endswith("boundary") for d in annot.description]]) == [2.0, 2.0]


def test_nihon_kohden_all_blocks(files):
    raw = eeg2edf.read_raw(files["nk"], block="all")
    assert raw.n_times == 500
    meta = get_sidecar(raw)
    assert [s["start"] for s in meta["segments"]] == ["2024-02-06T19:55:47",
                                                     "2024-02-06T19:56:49"]
    assert [s["offset_s"] for s in meta["segments"]] == [0.0, 2.0]
    annot = raw.annotations
    assert (3.0, "Seizure") in zip(annot.onset, annot.description)
    assert {d for d, o in zip(annot.description, annot.onset) if o == 2.0} >= {
        "BAD boundary", "EDGE boundary"}
    np.testing.assert_array_equal(
        raw.get_data()[:, 200:],
        eeg2edf.read_raw(files["nk"], block=1).get_data())


def test_fif_round_trip_keeps_the_sidecar(files, tmp_path):
    raw = eeg2edf.read_raw(files["nic"], preload=True)
    raw.save(tmp_path / "x_raw.fif", verbose="error")
    back = mne.io.read_raw_fif(tmp_path / "x_raw.fif", verbose="error")
    assert get_sidecar(back) == get_sidecar(raw)
    assert list(back.annotations.extras) == list(raw.annotations.extras)
    assert back.info["subject_info"] == raw.info["subject_info"]
    np.testing.assert_allclose(back.get_data(), raw.get_data(), rtol=1e-6)


def test_get_sidecar_refuses_other_raws():
    info = mne.create_info(["a"], 100.0)
    raw = mne.io.RawArray(np.zeros((1, 10)), info, verbose="error")
    with pytest.raises(ValueError, match="eeg2edf sidecar"):
        get_sidecar(raw)
    info["description"] = json.dumps({"schema": "something-else"})
    with pytest.raises(ValueError):
        get_sidecar(info)


# --- montages, channel types, positions ------------------------------------

def test_nihon_kohden_montage_matches_nk2edf(files, tmp_path):
    converted(nk2edf, files["nk"], tmp_path, "--blocks", "0", "--montage", "EMU1")
    ref = edf(next(tmp_path.glob("*.edf")))
    raw = eeg2edf.read_raw(files["nk"], montage="auto")
    assert raw.ch_names == ["Fp1-Fp2", "Fp2-C3", "0V-EKG1", "DC01", "Events/Markers"]
    meta = get_sidecar(raw)
    assert meta["montage_applied"] == "EMU1"
    # as nk2edf's sidecar has it: 0V-EKG1 derives from EKG1 against 0V
    assert [c["derived"] for c in meta["channels"][:3]] == [True, True, True]
    assert meta["channels"][2]["reference"] == "EKG1"
    traces = raw.get_data(picks=[0, 1, 2])
    # the EDF halves a bipolar trace into int16; one of its LSBs is the bound
    edf_lsb = np.array([c["resolution_uv_per_lsb"] for c in json.loads(
        next(tmp_path.glob("*.json")).read_text())["channels"][:3]])[:, None] * 1e-6
    assert_within(traces, ref.get_data(picks=[0, 1, 2]), edf_lsb)
    referential = eeg2edf.read_raw(files["nk"])
    np.testing.assert_allclose(traces[2], -referential.get_data(picks="EKG1")[0])
    assert raw.annotations.description.tolist() == \
        referential.annotations.description.tolist()


def test_micromed_montage(files):
    raw = eeg2edf.read_raw(files["vwr"], montage="Bipolar")
    assert raw.ch_names == ["A-B", "B-REF"]
    np.testing.assert_allclose(raw.get_data()[0], 0, atol=1e-15)
    np.testing.assert_allclose(raw.get_data()[1] * 1e6, [1, 2, 3, 4])
    # HISTORY recorded "Ground": A against its own ground
    assert eeg2edf.read_raw(files["vwr"], montage="auto").ch_names == ["A-REF"]


def test_nicolet_montage_after_fif(files, tmp_path):
    raw = eeg2edf.read_raw(files["nic"], preload=True)
    raw.save(tmp_path / "n_raw.fif", verbose="error")
    back = mne.io.read_raw_fif(tmp_path / "n_raw.fif", verbose="error")
    bip = apply_montage(back, "Bipolar")
    assert bip.ch_names == ["Fp1-Fp2", "EKG-Bipolar"]
    d = raw.get_data()
    np.testing.assert_allclose(bip.get_data()[0], d[0] - d[1], atol=1e-12)
    np.testing.assert_allclose(bip.get_data()[1], d[2], atol=1e-12)
    with pytest.raises(ValueError, match="no montage"):
        apply_montage(back, "nope")


def test_ch_types_and_positions(files):
    seeg = eeg2edf.read_raw(files["nk"], ch_types="seeg")
    assert seeg.get_channel_types() == ["seeg", "seeg", "seeg", "ecg", "misc", "stim"]
    typed = eeg2edf.read_raw(files["nk"], ch_types={"C3": "eog"})
    assert typed.get_channel_types()[2] == "eog"
    with pytest.raises(ValueError, match="not in the recording"):
        eeg2edf.read_raw(files["nk"], ch_types={"nope": "eeg"})
    from eeg2edf.mne._base import standard_montage_name

    placed = eeg2edf.read_raw(files["nk"], positions=standard_montage_name("1020"))
    assert not np.isnan(placed.info["chs"][0]["loc"][:3]).any()


def test_read_raw_rejects_unknown_extension(tmp_path):
    with pytest.raises(ValueError, match="don't know"):
        eeg2edf.read_raw(tmp_path / "x.foo")


def test_bad_block_is_a_value_error(files):
    with pytest.raises(ValueError, match="out of range"):
        eeg2edf.read_raw(files["nk"], block=5)
    with pytest.raises(ValueError, match="no montage"):
        eeg2edf.read_raw(files["nk"], montage="nope")


# --- the eeg2edf command ---------------------------------------------------

def run(*args):
    return subprocess.run([sys.executable, "-m", "eeg2edf", *map(str, args)],
                          capture_output=True, text=True)


def test_cli_dispatches_by_extension(files):
    assert "datablocks" in run(files["nk"], "--list").stdout
    assert "segments" in run(files["nic"], "--list").stdout
    assert "Hz" in run(files["vwr"], "--list").stdout
    assert "usage: eeg2edf" in run("--help").stdout
    assert eeg2edf.__version__ in run("--version").stdout
    bad = run(files["dir"] / "x.foo")
    assert bad.returncode != 0 and "don't know" in bad.stderr
