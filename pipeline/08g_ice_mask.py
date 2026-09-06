"""Ice-sheet domain mask for one 100 km tile from BedMachine Greenland v6 (mask: 0 ocean, 1 ice-free land, 2 grounded ice,
3 floating ice, 4 non-Greenland land; 150 m, EPSG:3413; Morlighem et al. 2017/2022, mask derived from the GIMP ice mask).
Writes out/{TILE}_icemask_150m.npy (bool, ice = grounded or floating) + out/{TILE}_icemask_150m_meta.json (x/y centres of
row 0 / col 0, step 150). 08e upsamples it (nearest) to the 10 m grid and uses it as the domain for water pixels.
Run:  TILE=29_45 $(cat .python_env) scripts/08g_ice_mask.py"""
import os, json, numpy as np, netCDF4 as nc
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(ROOT, "out")
TILE = os.environ.get("TILE", "19_39"); rr, cc = (int(v) for v in TILE.split("_"))
x0 = (cc - 1) * 100_000 - 4_000_000; y0 = (rr - 1) * 100_000 - 4_000_000; x1, y1 = x0 + 100_000, y0 + 100_000
BM = os.environ.get("BEDMACHINE", "/Users/jrines/stanford_gp/dissertation/defense/cinematic/source_data/data/BedMachineGreenland-v6.nc")
d = nc.Dataset(BM); x = d["x"][:]; y = d["y"][:]
c0, c1 = np.searchsorted(x, x0 - 150), np.searchsorted(x, x1 + 150); r0, r1 = np.searchsorted(-y, -(y1 + 150)), np.searchsorted(-y, -(y0 - 150))
m = np.asarray(d["mask"][r0:r1, c0:c1]); ice = (m == 2) | (m == 3)
np.save(os.path.join(OUT, f"{TILE}_icemask_150m.npy"), ice)
json.dump(dict(tile=TILE, x_col0=float(x[c0]), y_row0=float(y[r0]), step=150, shape=list(ice.shape), classes={int(k): int(v) for k, v in zip(*np.unique(m, return_counts=True))}),
          open(os.path.join(OUT, f"{TILE}_icemask_150m_meta.json"), "w"), indent=1)
print(f"{TILE}: window {ice.shape}, ice fraction {ice.mean():.3f}, classes {dict(zip(*np.unique(m, return_counts=True)))}")
