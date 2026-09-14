"""Relocate WITCH 1.0.10's runtime pointer; never write into site-packages at run time."""
import importlib.util
from pathlib import Path

package = Path(importlib.util.find_spec("witch_msa").origin).parent
config = package / "configs.py"
source = config.read_text()
old = "homepath = os.path.dirname(__file__) + '/home.path'"
new = """user_config = os.path.expanduser('~/.witch_msa')
os.makedirs(user_config, exist_ok=True)
homepath = os.path.join(user_config, 'home.path')"""
if source.count(old) != 1:
    raise RuntimeError("Unexpected WITCH configs.py: refusing to apply an unverified patch")
config.write_text(source.replace(old, new))
