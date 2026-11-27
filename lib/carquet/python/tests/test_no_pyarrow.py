"""The core tests pass with PyArrow unimportable."""

from __future__ import annotations

import os
import subprocess
import sys
import textwrap


def test_core_without_pyarrow(tmp_path):
    here = os.path.dirname(__file__)
    code = textwrap.dedent(
        f"""
        import sys
        sys.modules["pyarrow"] = None      # any `import pyarrow` now fails
        import pytest
        rc = pytest.main(["-q", "-p", "no:cacheprovider", "--basetemp", {str(tmp_path / 'bt')!r},
                          {os.path.join(here, 'test_core.py')!r}])
        loaded = [m for m, v in sys.modules.items() if v is not None and m.split(".")[0] == "pyarrow"]
        assert not loaded, loaded
        sys.exit(int(rc))
        """
    )
    import carquet

    env = dict(os.environ)
    pkg_root = os.path.dirname(os.path.dirname(os.path.abspath(carquet.__file__)))
    env["PYTHONPATH"] = os.pathsep.join([pkg_root, here] + [p for p in [env.get("PYTHONPATH")] if p])
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env, cwd=here)
    assert out.returncode == 0, out.stdout[-4000:] + out.stderr[-2000:]
