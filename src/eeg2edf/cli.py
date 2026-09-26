"""`eeg2edf INPUT ...`: run the converter that matches INPUT's extension.

Every option after INPUT is passed through unchanged, so
`eeg2edf FILE.EEG OUT --montage auto` is `nk2edf FILE.EEG OUT --montage auto`.
"""
import os
import sys

TOOLS = {
    ".eeg": ("nk2edf", "eeg2edf.nihon_kohden.convert"),
    ".e": ("nicolet2edf", "eeg2edf.nicolet.convert"),
    ".vwr": ("vwr2edf", "eeg2edf.micromed.convert"),
}

USAGE = """usage: eeg2edf INPUT [OUTDIR] [options]

Converts a clinical EEG recording to EDF+ with a JSON sidecar, choosing the
converter from INPUT's extension:

  .EEG   Nihon Kohden        (same as nk2edf)
  .e     Nicolet / Nervus    (same as nicolet2edf)
  .vwr   Micromed VWR        (same as vwr2edf)

`eeg2edf INPUT --help` lists that converter's options.
"""


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help"):
        print(USAGE)
        return
    if argv[0] in ("-V", "--version"):
        from . import __version__

        print(f"eeg2edf {__version__}")
        return
    ext = os.path.splitext(argv[0])[1].lower()
    if ext not in TOOLS:
        raise SystemExit(f"eeg2edf: don't know how to read {ext or 'a file with no extension'!r} "
                         f"files -- expected one of {', '.join(sorted(TOOLS))}")
    _name, module = TOOLS[ext]
    import importlib

    importlib.import_module(module).main(argv)


if __name__ == "__main__":
    main()
