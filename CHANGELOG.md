# Changelog

## Unreleased

- **Fix: `vwr2edf` on 32-bit recordings kept only 1/65,537 of the source
  resolution.** The whole 2^32 integer range was spread over EDF's 2^16 steps.
  4-byte channels now use their calibration range `lmin..lmax`, widened to the
  values actually stored so nothing is clipped. Re-convert any 4-byte VWR EDFs.
  1- and 2-byte output is unchanged.
- `vwr2edf --list` prints the sample width and, per channel, the EDF resolution
  and whether it is exact.
- README: a "Precision: what is lossless" section on which inputs and options
  round-trip exactly, and how to read that from the sidecar.

## 0.1.0

First release as a package.

- Installable from PyPI: `nk2edf`, `nicolet2edf`, `vwr2edf`, `eeg2edf` (picks
  one by extension) and `eeg2edf-bipolar-mtg` are console commands. The scripts
  no longer run from a checkout (`python nk2edf/nk2edf.py` is now `nk2edf`).
- MNE-Python readers (`pip install "eeg2edf[mne]"`): `eeg2edf.read_raw` and
  `eeg2edf.mne.read_raw_nihon_kohden` / `read_raw_nicolet` / `read_raw_micromed`
  return a lazily-read `mne.io.Raw`, with the sidecar in `info["description"]`,
  events as annotations and the vendor montages available through
  `apply_montage`.
- **Fix: `nk2edf --montage` wrote bipolar traces at half their amplitude.** The
  traces were halved into int16 as intended, but the EDF physical range was not
  doubled to match, so every montage trace read back at 50% of its true µV
  (and the sidecar's `resolution_uv_per_lsb` was half the true value).
  Referential output was unaffected, as were `nicolet2edf` and `vwr2edf`.
  Re-convert any Nihon Kohden montage EDFs made before this release.
- Output is otherwise byte-identical to the scripts this release replaces.
