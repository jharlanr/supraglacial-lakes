"""Ice-sheet domain mask for one 100 km tile from BedMachine Greenland v6 (mask: 0 ocean, 1 ice-free land, 2 grounded ice,
3 floating ice, 4 non-Greenland land; 150 m, EPSG:3413; Morlighem et al. 2017/2022, mask derived from the GIMP ice mask).
Writes out/{TILE}_icemask_150m.npy (bool, ice = grounded or floating) + out/{TILE}_icemask_150m_meta.json (x/y centres of
row 0 / col 0, step 150). 08e upsamples it (nearest) to the 10 m grid and uses it as the domain for water pixels.

Also writes out/{TILE}_extnonice_150m.npy (bool, same grid): the non-ice that is connected to the OUTSIDE of the
ice sheet, as opposed to nunataks and other non-ice enclosed by ice.  Spec 13c.1 uses it to tell an ice-marginal
lake (its outline reaches the ice sheet's edge) from a supraglacial lake that merely abuts a nunatak.  Connectivity
is flood-filled from the border of a window padded by PAD_KM (default 10 km) beyond the tile, so the answer does
not depend on where the 100 km tile happens to be cut; the result is then cropped back to the tile window.
Run:  TILE=29_45 $(cat .python_env) scripts/08g_ice_mask.py"""
import os, json, numpy as np, netCDF4 as nc, scipy.ndimage as ndi
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(ROOT, "out")
TILE = os.environ.get("TILE", "19_39"); rr, cc = (int(v) for v in TILE.split("_"))
x0 = (cc - 1) * 100_000 - 4_000_000; y0 = (rr - 1) * 100_000 - 4_000_000; x1, y1 = x0 + 100_000, y0 + 100_000
PAD = float(os.environ.get("PAD_KM", "10")) * 1000
BM = os.environ.get("BEDMACHINE", "/Users/jrines/stanford_gp/dissertation/defense/cinematic/source_data/data/BedMachineGreenland-v6.nc")
d = nc.Dataset(BM); x = d["x"][:]; y = d["y"][:]
c0, c1 = np.searchsorted(x, x0 - 150), np.searchsorted(x, x1 + 150); r0, r1 = np.searchsorted(-y, -(y1 + 150)), np.searchsorted(-y, -(y0 - 150))
m = np.asarray(d["mask"][r0:r1, c0:c1]); ice = (m == 2) | (m == 3)
np.save(os.path.join(OUT, f"{TILE}_icemask_150m.npy"), ice)

# --- non-ice connected to the outside of the ice sheet (padded, so it is tile-cut independent) ---
pc0, pc1 = np.searchsorted(x, x0 - PAD), np.searchsorted(x, x1 + PAD)
pr0, pr1 = np.searchsorted(-y, -(y1 + PAD)), np.searchsorted(-y, -(y0 - PAD))
pm = np.asarray(d["mask"][pr0:pr1, pc0:pc1]); pice = (pm == 2) | (pm == 3)
plab, _ = ndi.label(~pice)
edge_ids = set(np.unique(np.concatenate([plab[0], plab[-1], plab[:, 0], plab[:, -1]]))); edge_ids.discard(0)
pext = np.isin(plab, list(edge_ids))
ext = pext[r0 - pr0:r0 - pr0 + ice.shape[0], c0 - pc0:c0 - pc0 + ice.shape[1]]   # crop back to the tile window
assert ext.shape == ice.shape, (ext.shape, ice.shape)
np.save(os.path.join(OUT, f"{TILE}_extnonice_150m.npy"), ext)
nonice = ~ice
print(f"{TILE}: non-ice {nonice.sum()*150**2*1e-6:.0f} km2 = exterior {(nonice & ext).sum()*150**2*1e-6:.0f} + "
      f"interior (nunataks) {(nonice & ~ext).sum()*150**2*1e-6:.0f}; pad {PAD/1000:.0f} km")
json.dump(dict(tile=TILE, x_col0=float(x[c0]), y_row0=float(y[r0]), step=150, shape=list(ice.shape), classes={int(k): int(v) for k, v in zip(*np.unique(m, return_counts=True))}),
          open(os.path.join(OUT, f"{TILE}_icemask_150m_meta.json"), "w"), indent=1)
print(f"{TILE}: window {ice.shape}, ice fraction {ice.mean():.3f}, classes {dict(zip(*np.unique(m, return_counts=True)))}")
