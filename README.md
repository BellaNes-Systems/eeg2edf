# eeg2edf

Converters from proprietary clinical EEG formats to **EDF+**, plus a JSON sidecar
carrying the metadata EDF+ has nowhere to put.

| Tool | Reads | Notes |
|---|---|---|
| `nk2edf/` | Nihon Kohden EEG-1100 / EEG-1200A (`.EEG` + `.21E`/`.PNT`/`.LOG`/`.PTN`) | one EDF+C per clip; `.LOG` events become annotations |
| `nicolet2edf/` | Nicolet / Nervus `.e` | one EDF+C per stored segment; `--concat` for the vendor timeline |
| `vwr2edf/` | Micromed VWR / Brain-Quick `.vwr` | montage, events and trigger tracks decoded |

All three are lossless: the EDF digital range is the stored range, so sample
values round-trip exactly.

## Why

These formats are readable only through vendor software that is discontinued,
Windows-only, or dependent on COM DLLs that no longer ship. Each converter here
was written by reverse-engineering the container against real recordings and
verifying the output independently — see each tool's README for the evidence.

## Install

```bash
pip install numpy          # the only runtime dependency
pip install edfio          # optional, for round-trip validation in tests
```

No packaging step: run the scripts in place.

## Use

```bash
python nk2edf/nk2edf.py       INPUT.EEG OUTDIR
python nicolet2edf/nicolet2edf.py INPUT.e   OUTDIR
python vwr2edf/vwr2edf.py     INPUT.vwr OUTDIR
```

Every tool takes `--list` to report what a file contains without converting, and
`--no-sidecar` to skip the JSON. Per-tool options are in the tool's README.

`make_bipolar_mtg.py` is a small extra: it emits an EDFbrowser `.mtg` montage
from an SEEG EDF.

## The sidecar

Each conversion writes `OUTPUT.json` beside `OUTPUT.edf` following the
`eeg2edf-sidecar/1` schema — clip, patient, channel, trace, event and segment
records. **[SIDECAR.md](SIDECAR.md) is the normative spec.** `edfcommon.py`
builds it and is shared by all three converters.

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
python test_sidecar_schema.py        # all three converters against the schema
cd vwr2edf && python test_vwr2edf.py
```

## Licensing

Apache-2.0 (`LICENSE`), **except `vwr2edf/`**, which is BSD-3-Clause
(`vwr2edf/LICENSE`) because it contains portions derived from libvwr,
Copyright (C) Franco Milicchio. That notice must be retained when
redistributing those portions. See `NOTICE`.
