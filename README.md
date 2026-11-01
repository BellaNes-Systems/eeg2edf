# eeg2edf

Readers for proprietary clinical EEG formats: converters to **EDF+**
with a JSON sidecar carrying the metadata EDF+ has nowhere to put, and
**MNE-Python** readers that open the same recordings directly as `mne.io.Raw`.

| Format | Command | MNE reader | Notes |
|---|---|---|---|
| Nihon Kohden EEG-1100 / EEG-1200A (`.EEG` + `.21E`/`.PNT`/`.LOG`/`.PTN`) | `nk2edf` | `read_raw_nihon_kohden` | one EDF+C per clip; `.LOG` events become annotations |
| Nicolet / Nervus `.e` | `nicolet2edf` | `read_raw_nicolet` | one EDF+C per stored segment; `--concat` for the vendor timeline |
| Micromed VWR / Brain-Quick `.vwr` | `vwr2edf` | `read_raw_micromed` | montage, events and trigger tracks decoded |

Referential output from 8- and 16-bit sources is lossless: the EDF digital
range is the stored range, so sample values round-trip exactly. 32-bit Micromed
samples and `--montage` traces have to be re-quantised to fit EDF's 16 bits —
see [Precision](#precision-what-is-lossless). The MNE readers are exact in every
case: they hand MNE the stored integers with their exact calibration.

## Why

These formats are readable only through vendor software that is discontinued,
Windows-only, or dependent on COM DLLs that no longer ship. Each reader here
was written by reverse-engineering the container against real recordings and
verifying the output independently — see each tool's notes in [docs/](https://github.com/BellaNes-Systems/eeg2edf/tree/main/docs)
for the evidence.

## Install

```bash
pip install eeg2edf            # converters; numpy is the only dependency
pip install "eeg2edf[mne]"     # plus the MNE-Python readers
```

## Convert to EDF+

```bash
nk2edf      INPUT.EEG OUTDIR
nicolet2edf INPUT.e   OUTDIR
vwr2edf     INPUT.vwr OUTDIR
eeg2edf     INPUT     OUTDIR      # picks the converter from the extension
```

Every tool takes `--list` to report what a file contains without converting, and
`--no-sidecar` to skip the JSON. Per-tool options are in
[docs/nk2edf.md](https://github.com/BellaNes-Systems/eeg2edf/blob/main/docs/nk2edf.md), [docs/nicolet2edf.md](https://github.com/BellaNes-Systems/eeg2edf/blob/main/docs/nicolet2edf.md) and
[docs/vwr2edf.md](https://github.com/BellaNes-Systems/eeg2edf/blob/main/docs/vwr2edf.md).

`eeg2edf-bipolar-mtg` is a small extra: it emits an EDFbrowser `.mtg` montage
from an SEEG EDF.

## Precision: what is lossless

EDF stores every sample as a 16-bit integer (an LSB is one integer step). An
output is exact when each stored value gets its own EDF step. It is quantised
when the source needs more than 16 bits, and then the EDF keeps the nearest
step, within half an output LSB.

| Input / option | EDF output |
|---|---|
| Nihon Kohden, referential (EEG, DC, event word) | exact |
| Nicolet, referential | exact. The default polarity flip maps a stored −32768, which has no positive int16, to +32767; `--no-invert` is exact for every value |
| Micromed, 8- or 16-bit samples, referential | exact |
| Micromed, 32-bit samples | quantised: the channel's calibration range `lmin..lmax`, widened to the values actually stored, spread over 65,536 steps. Nothing is clipped |
| `--montage`, a trace that subtracts two channels (all three tools) | quantised: the difference needs 17 bits, so the output LSB is twice the source LSB |
| `--montage`, a single-channel trace | Nicolet (`EKG-Bipolar`) and Micromed (ground-referenced, `ECG1+`): exact. **Nihon Kohden: halved like the rest** (`-EKG1`) |
| MNE readers, `apply_montage` | exact; montages are computed in float64 and never re-quantised |

The sidecar shows the boundary per channel. `resolution_uv_per_lsb` is the EDF
step in µV, in full precision; the EDF header holds the same gain in 8 ASCII
characters, rounded to within a few percent of an LSB. A channel is exact when
that step is no larger than the source's own: the `factor` field for Micromed,
the referential resolution for Nihon Kohden and Nicolet, where `derived: true`
traces are twice as coarse. To see it before converting, `nicolet2edf --list`
prints each channel's resolution, and `vwr2edf --list` prints the sample width
and, per channel, the EDF step and `exact` or how many times coarser it is.

## Read with MNE-Python

```python
import eeg2edf
from eeg2edf.mne import apply_montage, get_sidecar, montage_names

raw = eeg2edf.read_raw("FILE.EEG")                # .EEG, .e or .vwr; data read lazily
raw = eeg2edf.read_raw("FILE.EEG", block="all")   # every Nihon Kohden clip, joined
raw = eeg2edf.read_raw("SAMPLE.e", segment=2)     # one Nicolet segment

raw.annotations                   # every event, with its source fields as .extras
meta = get_sidecar(raw)           # the eeg2edf-sidecar/1 dict the converter writes

montage_names(raw)                # the vendor display montages, e.g. ['ETEST', 'EMU1']
bipolar = apply_montage(raw, "EMU1")
raw = eeg2edf.read_raw("FILE.EEG", montage="auto")  # the montage the .LOG names
```

What MNE has a field for goes there:

| From the recording | In MNE |
|---|---|
| clip start | `info["meas_date"]` (vendor wall clock, labelled UTC as MNE's EDF reader does) |
| sex, birth date | `info["subject_info"]` — never a name or record number |
| device | `info["device_info"]` |
| filter settings, notch (Nicolet) | `info["highpass"]`, `info["lowpass"]`, `info["line_freq"]` |
| channel units and resolution | each channel's `cal` and unit; ECG/EOG/EMG, DC (misc) and the NK event word (stim) typed |
| events | `raw.annotations`, with type, source, GUID and timestamp in `extras` and the event's channel in `ch_names` |
| segment joins (Nicolet, NK `block="all"`, VWR acquisition cuts) | `BAD boundary` / `EDGE boundary` annotations |

Everything else — the vendor montages, per-channel references and filters,
segment start times, age, events that fall between segments — is the sidecar
itself, stored as JSON in `info["description"]`. It survives `copy()`, `crop()`,
`pick()` and a FIF save/load, and `get_sidecar(raw)` reads it back. It
describes the channels as they were read: picking or renaming channels later
does not rewrite it.

Other reader options: `ch_types="seeg"` (or a dict) overrides the guessed
channel types; `positions="standard_1020"` (`"colin27_1020"` from MNE 1.13)
sets sensor positions — off by default, because SEEG contact names such as `A1`
collide with 10-20 names; `include_trends=True` adds Nicolet's derived trend
channels, held onto the EEG rate. See each reader's docstring.

## The sidecar

Each conversion writes `OUTPUT.json` beside `OUTPUT.edf` following the
`eeg2edf-sidecar/1` schema — clip, patient, channel, trace, event and segment
records. **[SIDECAR.md](https://github.com/BellaNes-Systems/eeg2edf/blob/main/SIDECAR.md) is the normative spec.**
`eeg2edf.edfcommon` builds it and is shared by all three converters and the
MNE readers.

The schema has an external consumer:
[bellanes-lab](https://github.com/BellaNes-Systems/bellanes-lab) reads it to
recover per-file montage and reference information. Nothing imports across the
two repositories — the contract is the JSON — but changing the schema means
checking that reader.

The sidecar never carries a patient name or medical record number. The EDF+
patient field is limited to sex and birth date, which affect interpretation;
`--patient` overrides it.

## Tests

Fully synthetic — no recordings needed.

```bash
pip install -e ".[test]"
pytest
```

`tests/synth.py` writes small recordings in all three formats from the layouts
the readers document; the MNE tests hold each reader to its converter's EDF.

## Releasing

Releases go to PyPI from GitHub Actions (`.github/workflows/release.yml`) by
[trusted publishing](https://docs.pypi.org/trusted-publishers/) — no API token
is stored anywhere.

One-time setup:

1. On [PyPI](https://pypi.org/manage/account/publishing/) and
   [TestPyPI](https://test.pypi.org/manage/account/publishing/), add a pending
   publisher: project `eeg2edf`, owner `BellaNes-Systems`, repository `eeg2edf`,
   workflow `release.yml`, environment `pypi` (TestPyPI: `testpypi`).
2. In the GitHub repository settings, create the environments `pypi` and
   `testpypi` (optionally requiring a reviewer for `pypi`).

Each release: bump `__version__` in `src/eeg2edf/__init__.py`, add a
[CHANGELOG.md](https://github.com/BellaNes-Systems/eeg2edf/blob/main/CHANGELOG.md) entry, merge, then publish a GitHub Release tagged
`v<version>`. The workflow checks the tag matches, builds, uploads to TestPyPI,
then to PyPI.

## Licensing

Apache-2.0 (`LICENSE`), **except `src/eeg2edf/micromed/`**, which is
BSD-3-Clause (`src/eeg2edf/micromed/LICENSE`) because it contains portions
derived from libvwr, Copyright (C) Franco Milicchio. That notice must be
retained when redistributing those portions, and ships in the wheel. See
`NOTICE`.
