import sys
from pathlib import Path

# synth.py is a plain helper module next to the tests, not part of the package.
sys.path.insert(0, str(Path(__file__).parent))
