"""Shared paths, seeds and helpers for paper_b scripts. Read-only w.r.t. core-data sources."""
import json
import platform
import random
import sys
from datetime import datetime
from importlib import metadata
from pathlib import Path

import numpy as np

PAPER = Path(__file__).resolve().parent
CORE = PAPER.parent
RHAI = CORE / "RHAI"
SWAI = CORE / "SWAI"
FLOOD = CORE / "flood-recovery-sindh"
PAPER_A = CORE / "paper_a"
DATA = PAPER / "data"
TABLES = PAPER / "results" / "tables"
FIGS = PAPER / "results" / "figures"
DOCS = PAPER / "docs"
LOGS = PAPER / "logs"
for _d in (DATA, TABLES, FIGS, LOGS):
    _d.mkdir(parents=True, exist_ok=True)

SEED = 20221001
VID = "village_id"
CID = "cell_id"
PKGS = ["numpy", "pandas", "pyarrow", "scipy", "statsmodels", "linearmodels", "scikit-learn", "econml",
        "geopandas", "shapely", "pyproj", "pyogrio", "rasterio", "networkx", "osmnx", "wildboottest", "matplotlib"]


def set_seeds(seed=SEED):
    random.seed(seed)
    np.random.seed(seed)


def log_env(script_name):
    """Write package versions + seed for this run to logs/<script>.env.json."""
    vers = {}
    for p in PKGS:
        try:
            vers[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            vers[p] = None
    info = {"script": Path(script_name).name, "run_at": datetime.now().isoformat(timespec="seconds"),
            "python": sys.version.split()[0], "platform": platform.platform(),
            "seed": SEED, "packages": vers}
    (LOGS / f"{Path(script_name).stem}.env.json").write_text(json.dumps(info, indent=2))
    return info
