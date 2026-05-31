"""Make the ``src`` layout importable during local development without an install.

A normal ``pip install -e .`` puts ``src`` on the path via the package's editable record.
This shim is a belt-and-braces fallback so the test suite (and ``python conftest.py``-style
ad-hoc runs) work even when the editable install isn't present. It is intentionally inert
when ``openth`` is already importable.
"""

import sys
from pathlib import Path

_SRC = Path(__file__).parent / "src"
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))
