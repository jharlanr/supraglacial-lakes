"""Site-level persistence test: build per-year outlines from the 08b count rasters for every season on
disk (TOA NDWI>0.5, k=1, closing 5 px, no opening), then sites = union over years with a closing, and score
against Dunmire's lakes of every year. Variants:
  V1  per-year minimum 0.05 km2 + fill>=0.5, site = union of per-year components
  V2  per-year minimum 0.01 km2 + fill>=0.5, site = union, keep sites with area >= 0.05 km2 AND seen in >= 2 years
  V3  as V2 but seen in >= 1 year (isolates the effect of the recurrence requirement)
Writes out/08_multiyear_sites_{TILE}.{txt,csv,png}.
Run:  nice -n 15 $(cat .python_env) scripts/08d_multiyear_sites.py        (TILE=19_39; SITE_CLOSE_PX=15)
"""
import os, json, glob, time
import numpy as np, pandas as pd, geopandas as gpd
from scipy import ndimage as ndi
from rasterio.features import rasterize
from rasterio.transform import Affine
from shapely.geometry import box
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(ROOT, "out")
TILE = os.environ.get("TILE", "19_39"); SITE_CLOSE = int(os.environ.get("SITE_CLOSE_PX", "15")); FILL = 0.5
years = sorted(int(os.path.basename(f)[len(f"08_s2counts_{TILE}_"):].split("_")[0]) for f in glob.glob(os.path.join(OUT, f"08_s2counts_{TILE}_*_meta.json")))
meta = json.load(open(os.path.join(OUT, f"08_s2counts_{TILE}_{years[0]}_meta.json"))); tr = Affine(*meta["transform"]); N = meta["shape"][0]
print("seasons on disk:", years)
def disk(r):
    y, x = np.ogrid[-r:r + 1, -r:r + 1]; return (x * x + y * y) <= r * r
S8 = np.ones((3, 3), bool)

def year_components(arr, min_px):
    w = ndi.binary_closing(arr >= 1, structure=disk(5)); lab, _ = ndi.label(w, structure=S8)
    cnt = np.bincount(lab.ravel()); ids = np.flatnonzero(cnt >= min_px); ids = ids[ids > 0]
    fill = ndi.mean((arr >= 1).astype(float), lab, ids); keep = ids[fill >= FILL]
    return np.where(np.isin(lab, keep), lab, 0)

counts = {y: np.load(os.path.join(OUT, f"08_s2counts_{TILE}_{y}_n_w50_toa.npy")) for y in years}
dun = {}; dlab = {}; dpx = {}
tile_box = box(meta["x0"], meta["y0"], meta["x1"], meta["y1"])
for y in years:
    g = gpd.read_file(os.path.join(ROOT, "..", "labels", "dunmire", f"labels_{y}_volumes.geojson")).to_crs(3413)
    dun[y] = g[g.geometry.within(tile_box)].reset_index(drop=True)
    dlab[y] = rasterize(((geom, i + 1) for i, geom in enumerate(dun[y].geometry)), out_shape=(N, N), transform=tr, dtype="int32")
    dpx[y] = np.bincount(dlab[y].ravel())[1:]
    print(f"Dunmire {y}: {len(dun[y])} lakes inside the tile")

def score(site_lab, name):
    sites = np.unique(site_lab[site_lab > 0]); spx = np.bincount(site_lab.ravel())
    out = dict(variant=name, n_sites=len(sites), site_km2=round(spx[sites].sum() * 1e-4, 1))
    touched_any = set()
    for y in years:
        m = site_lab > 0
        cov = np.bincount(dlab[y][(dlab[y] > 0) & m], minlength=len(dun[y]) + 1)[1:] / np.maximum(dpx[y], 1)
        out[f"dun{y}_touched"] = int((cov > 0).sum()); out[f"dun{y}_matched"] = int((cov >= 0.5).sum()); out[f"dun{y}_total"] = len(dun[y])
        touched_any |= set(np.unique(site_lab[(site_lab > 0) & (dlab[y] > 0)]))
    extras = np.setdiff1d(sites, list(touched_any))
    out["extras"] = len(extras); out["extra_km2"] = round(spx[extras].sum() * 1e-4, 1)
    # lid test: Dunmire lakes of different years that share one site with zero polygon overlap is not computed here; report sites holding >=2 lakes of one year
    for y in years:
        pairs = np.unique(np.stack([site_lab[(site_lab > 0) & (dlab[y] > 0)], dlab[y][(site_lab > 0) & (dlab[y] > 0)]], 1), axis=0)
        per = pd.Series(pairs[:, 1]).groupby(pairs[:, 0]).nunique() if len(pairs) else pd.Series(dtype=int)
        out[f"sites_with_ge2_lakes_{y}"] = int((per >= 2).sum()); out[f"max_lakes_per_site_{y}"] = int(per.max()) if len(per) else 0
    return out

rows = []; keep_labs = {}
for name, min_px, min_years in (("V1 per-year 0.05 km2, union", 500, 1), ("V2 per-year 0.01 km2, union >=0.05 km2, >=2 yrs", 100, 2), ("V3 per-year 0.01 km2, union >=0.05 km2, >=1 yr", 100, 1)):
    ylabs = {y: year_components(counts[y], min_px) for y in years}
    union = np.zeros((N, N), bool); nyears = np.zeros((N, N), np.uint8)
    for y in years: union |= ylabs[y] > 0
    u = ndi.binary_closing(union, structure=disk(SITE_CLOSE)) if SITE_CLOSE else union
    slab, _ = ndi.label(u, structure=S8); spx = np.bincount(slab.ravel()); ids = np.flatnonzero(spx >= 500); ids = ids[ids > 0]
    # years in which each site holds a per-year component
    seen = np.zeros(slab.max() + 1, int)
    for y in years:
        present = np.unique(slab[(slab > 0) & (ylabs[y] > 0)]); seen[present] += 1
    keep = ids[seen[ids] >= min_years]; slab = np.where(np.isin(slab, keep), slab, 0)
    r = score(slab, name); rows.append(r); keep_labs[name] = slab; print(r)
df = pd.DataFrame(rows); df.to_csv(os.path.join(OUT, f"08_multiyear_sites_{TILE}.csv"), index=False)

# which of the 2019 misses (untouched at k=1 single-year) are recovered by each variant
base19 = year_components(counts[2019], 500) if 2019 in counts else None
lines = [f"Multi-year sites, tile {TILE}, seasons {years}, site closing {SITE_CLOSE} px — {time.strftime('%Y-%m-%d %H:%M %Z')}", df.to_string(index=False)]
if base19 is not None:
    cov = np.bincount(dlab[2019][(dlab[2019] > 0) & (base19 > 0)], minlength=len(dun[2019]) + 1)[1:] / np.maximum(dpx[2019], 1)
    missed = np.flatnonzero(cov == 0); lines.append(f"\n2019 Dunmire lakes untouched by the single-season recipe: {len(missed)}")
    for name, slab in keep_labs.items():
        c2 = np.bincount(dlab[2019][(dlab[2019] > 0) & (slab > 0)], minlength=len(dun[2019]) + 1)[1:] / np.maximum(dpx[2019], 1)
        lines.append(f"  recovered by {name}: {int((c2[missed] > 0).sum())} touched, {int((c2[missed] >= 0.5).sum())} covered >= 50 %")
open(os.path.join(OUT, f"08_multiyear_sites_{TILE}.txt"), "w").write("\n".join(lines) + "\n"); print("\n".join(lines[2:]))

# figure: crescent, dumbbell, wide view — V2 sites with Dunmire 2018 (orange) and 2019 (navy)
views = {"crescent CW_0381": (-172_900, -2_144_700, 1_500), "dumbbell CW2018_1270": (-185_200, -2_109_500, 2_200), "wide view": (meta["x0"] + 50_000, meta["y0"] + 50_000, 12_000)}
fig, axes = plt.subplots(len(views), len(keep_labs), figsize=(5.6 * len(keep_labs), 5.4 * len(views)))
for i, (vn, (cx, cy, half)) in enumerate(views.items()):
    c0, r0 = ~tr * (cx - half, cy + half); c1, r1 = ~tr * (cx + half, cy - half); r0, r1, c0, c1 = int(max(r0, 0)), int(min(r1, N)), int(max(c0, 0)), int(min(c1, N)); ext = (cx - half, cx + half, cy - half, cy + half)
    for j, (name, slab) in enumerate(keep_labs.items()):
        ax = axes[i, j]; sub = slab[r0:r1, c0:c1]
        ax.imshow(np.where(sub > 0, (sub % 7) + 1, 0), extent=ext, cmap="Pastel1", vmin=0, vmax=8, interpolation="nearest")
        for y, col in ((2018, "tab:orange"), (2019, "navy")):
            if y in dun:
                for g in dun[y][dun[y].geometry.intersects(box(ext[0], ext[2], ext[1], ext[3]))].geometry:
                    xx, yy = g.exterior.xy; ax.plot(xx, yy, color=col, lw=1)
        ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(f"{vn}: {name}", fontsize=8)
fig.suptitle(f"Sites from the union of {years} (pastel = site); orange = Dunmire 2018, navy = 2019", fontsize=10); fig.tight_layout()
fig.savefig(os.path.join(OUT, f"08_multiyear_sites_{TILE}.png"), dpi=95); print("wrote", f"out/08_multiyear_sites_{TILE}.{{txt,csv,png}}")
