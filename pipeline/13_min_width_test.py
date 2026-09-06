"""Minimum-width test for per-season outlines (Josh, 2026-09-06): open each outline with a disk of R px (removes anything
narrower than 2R px and cuts pond chains at necks), then grow the survivors back by R px inside the original outline
(shapes restored, channels gone). Reports outlines, area, and Dunmire recall per season with and without the rule, and
draws a before/after window.   Run:  TILE=29_45 nice -n 15 $(cat .python_env) scripts/13_min_width_test.py  (R_PX=3 YEARS="2018 2019")"""
import os, json, glob, numpy as np, pandas as pd, geopandas as gpd, scipy.ndimage as ndi, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from rasterio.features import rasterize; from affine import Affine; from shapely.geometry import box
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(ROOT, "out"); TILE = os.environ.get("TILE", "29_45"); R = int(os.environ.get("R_PX", "3"))
YEARS = [int(y) for y in os.environ.get("YEARS", "2018 2019").split()]; MIN_PX = 500; FILL = 0.5; S8 = np.ones((3, 3), bool)
def disk(r):
    y, x = np.ogrid[-r:r + 1, -r:r + 1]; return (x * x + y * y) <= r * r
pref = f"08_s2counts_{TILE}_"; meta = json.load(open(os.path.join(OUT, f"{pref}{YEARS[0]}_meta.json"))); tr = Affine(*meta["transform"]); N = meta["shape"][0]; x0, y0, x1, y1 = meta["x0"], meta["y0"], meta["x1"], meta["y1"]
im = np.load(os.path.join(OUT, f"{TILE}_icemask_150m.npy")); imm = json.load(open(os.path.join(OUT, f"{TILE}_icemask_150m_meta.json")))
xc = x0 + 10 * (np.arange(N) + 0.5); yc = y1 - 10 * (np.arange(N) + 0.5)
ice = im[np.clip(np.rint((imm["y_row0"] - yc) / 150).astype(int), 0, im.shape[0] - 1)[:, None], np.clip(np.rint((xc - imm["x_col0"]) / 150).astype(int), 0, im.shape[1] - 1)[None, :]]
dun = {y: gpd.read_file(os.path.join(ROOT, "..", "labels", "dunmire", f"labels_{y}_volumes.geojson")).to_crs(3413) for y in (2018, 2019)}
dun = {y: d[d.geometry.within(box(x0, y0, x1, y1))] for y, d in dun.items()}
def outlines(arr):
    w = ndi.binary_closing(arr >= 1, structure=disk(5)); lab, _ = ndi.label(w, structure=S8); cnt = np.bincount(lab.ravel()); ids = np.flatnonzero(cnt >= MIN_PX); ids = ids[ids > 0]
    keep = ids[ndi.mean((arr >= 1).astype(float), lab, ids) >= FILL]; keep = keep[ndi.maximum(ndi.binary_erosion(np.isin(lab, keep), structure=disk(2)).astype(np.uint8), lab, keep).astype(bool)]
    return np.isin(lab, keep)
def min_width(M):
    o = ndi.binary_opening(M, structure=disk(R)); g = ndi.binary_dilation(o, structure=disk(R)) & M
    lab, _ = ndi.label(g, structure=S8); cnt = np.bincount(lab.ravel()); ids = np.flatnonzero(cnt >= MIN_PX); ids = ids[ids > 0]; return np.isin(lab, ids)
rows = []; keepM = {}
for y in YEARS:
    arr = np.where(ice, np.load(os.path.join(OUT, f"{pref}{y}_n_w50_toa.npy")), 0); M0 = outlines(arr); M1 = min_width(M0); keepM[y] = (M0, M1)
    for name, M in (("current", M0), (f"min width {2*R*10} m", M1)):
        lab, n = ndi.label(M, structure=S8); row = dict(year=y, rule=name, outlines=n, area_km2=round(M.sum() * 1e-4, 1))
        if y in dun and len(dun[y]):
            dl = rasterize(((g, i + 1) for i, g in enumerate(dun[y].geometry)), out_shape=(N, N), transform=tr, dtype="int32"); hit = np.unique(dl[M]); hit = hit[hit > 0]
            row["dunmire_touched"] = f"{len(hit)}/{len(dun[y])}"
        rows.append(row)
t = pd.DataFrame(rows); print(t.to_string(index=False)); t.to_csv(os.path.join(OUT, f"13_min_width_{TILE}.csv"), index=False)
# before/after figure on the window with the most removed water
y = YEARS[-1]; M0, M1 = keepM[y]; diff = M0 & ~M1; k = 1200  # 12 km blocks
best = max(((i, j, diff[i:i + k, j:j + k].sum()) for i in range(0, N - k, 600) for j in range(0, N - k, 600)), key=lambda t: t[2])
i, j = best[:2]
if os.environ.get("CENTER"):  # "x,y" in EPSG:3413: centre the window there instead
    cx_, cy_ = (float(v) for v in os.environ["CENTER"].split(",")); c, r = ~tr * (cx_, cy_); i, j = int(max(min(r - k / 2, N - k), 0)), int(max(min(c - k / 2, N - k), 0))
ext = ((tr * (j, 0))[0], (tr * (j + k, 0))[0], (tr * (0, i + k))[1], (tr * (0, i))[1])
fig, axes = plt.subplots(1, 2, figsize=(14, 7.3))
for ax, M, ttl in zip(axes, (M0, M1), ("current outlines", f"with minimum width {2*R*10} m")):
    ax.imshow(M[i:i + k, j:j + k], extent=ext, cmap="Blues", vmin=0, vmax=1.4); ax.set_title(f"{TILE} {y}: {ttl}"); ax.set_xticks([]); ax.set_yticks([])
axes[1].contour(np.flipud(diff[i:i + k, j:j + k]), levels=[0.5], colors=["red"], linewidths=0.6, extent=ext, origin="lower"); axes[1].text(0.02, 0.02, "red = removed", color="red", transform=axes[1].transAxes)
fig.tight_layout(); fig.savefig(os.path.join(OUT, f"13_min_width_{TILE}{os.environ.get('TAG', '')}.png"), dpi=90); print("figure window rows", i, "cols", j, "removed px", best[2])
