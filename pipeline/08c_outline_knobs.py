"""Sweep the outline knobs on the downloaded 10 m water-detection counts (08b) and compare with
Dunmire's 2019 outlines inside the tile. Knobs: persistence k (scenes with NDWI>0.5), river opening
(0 or 3 px), closing radius (px), threshold band (n_w50_toa | n_w30_toa | n_w50_sr).
Writes out/08_outline_knobs_{TILE}_{YEAR}.{txt,csv,png}.

Run:  nice -n 15 $(cat .python_env) scripts/08c_outline_knobs.py       (TILE=19_39 YEAR=2019)
"""
import os, json, time, itertools
import numpy as np, pandas as pd, geopandas as gpd
from scipy import ndimage as ndi
from rasterio.features import rasterize, shapes
from rasterio.transform import Affine
from shapely.geometry import shape, box
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(ROOT, "out")
TILE = os.environ.get("TILE", "19_39"); YEAR = int(os.environ.get("YEAR", "2019"))
MIN_PX = 500                         # 0.05 km2 at 10 m
meta = json.load(open(os.path.join(OUT, f"08_s2counts_{TILE}_{YEAR}_meta.json")))
tr = Affine(*meta["transform"]); N = meta["shape"][0]
ld = lambda b: np.load(os.path.join(OUT, f"08_s2counts_{TILE}_{YEAR}_{b}.npy"))
n_obs = ld("n_obs"); bands = {b: ld(b) for b in ("n_w50_toa", "n_w30_toa", "n_w50_sr")}

# Dunmire 2019 outlines inside the tile, rasterised and labelled on the same grid
tile_box = box(meta["x0"], meta["y0"], meta["x1"], meta["y1"])
dun = gpd.read_file(os.path.join(ROOT, "..", "labels", "dunmire", f"labels_{YEAR}_volumes.geojson")).to_crs(3413)
dun = dun[dun.geometry.within(tile_box)].reset_index(drop=True)
dun_lab = rasterize(((g, i + 1) for i, g in enumerate(dun.geometry)), out_shape=(N, N), transform=tr, dtype="int32")
dun_px = np.bincount(dun_lab.ravel())[1:]
print(f"Dunmire {YEAR} lakes fully inside tile {TILE}: {len(dun)}, total {dun_px.sum()*1e-4:.1f} km2")

def disk(r):
    y, x = np.ogrid[-r:r + 1, -r:r + 1]; return (x * x + y * y) <= r * r

def run(band, k, open_px, close_r):
    w = bands[band] >= k
    if open_px: w = ndi.binary_opening(w, structure=np.ones((open_px, open_px), bool))
    if close_r: w = ndi.binary_closing(w, structure=disk(close_r))
    lab, n = ndi.label(w, structure=np.ones((3, 3), bool))
    cnt = np.bincount(lab.ravel()); keep = np.flatnonzero(cnt >= MIN_PX); keep = keep[keep > 0]
    m = np.isin(lab, keep); lab = np.where(m, lab, 0)
    # per Dunmire lake: fraction covered by any kept component; matched if >= 0.5
    cov = np.bincount(dun_lab[(dun_lab > 0) & m], minlength=len(dun) + 1)[1:] / np.maximum(dun_px, 1)
    matched = int((cov >= 0.5).sum()); touched = int((cov > 0).sum())
    # components: which touch a Dunmire lake; chaining = max number of Dunmire lakes under one component
    comp_ids = np.unique(lab[lab > 0]); comp_px = np.bincount(lab.ravel())
    pairs = np.unique(np.stack([lab[(lab > 0) & (dun_lab > 0)], dun_lab[(lab > 0) & (dun_lab > 0)]], 1), axis=0)
    per_comp = pd.Series(pairs[:, 1]).groupby(pairs[:, 0]).nunique() if len(pairs) else pd.Series(dtype=int)
    extra = np.setdiff1d(comp_ids, pairs[:, 0] if len(pairs) else [])
    return dict(band=band, k=k, open_px=open_px, close_r=close_r, n_comp=len(comp_ids),
                area_km2=round(comp_px[comp_ids].sum() * 1e-4, 1), dun_matched=matched, dun_touched=touched,
                dun_total=len(dun), n_extra=len(extra), extra_km2=round(comp_px[extra].sum() * 1e-4, 1),
                extra_median_km2=round(float(np.median(comp_px[extra])) * 1e-4, 3) if len(extra) else 0,
                max_chain=int(per_comp.max()) if len(per_comp) else 0,
                n_chain_ge3=int((per_comp >= 3).sum()) if len(per_comp) else 0,
                largest_km2=round(comp_px[comp_ids].max() * 1e-4, 2) if len(comp_ids) else 0), lab

rows = []; labs = {}
for band, k, op, cr in itertools.product(["n_w50_toa"], [1, 2, 3, 5], [0, 3], [5]):
    r, lab = run(band, k, op, cr); rows.append(r); labs[(band, k, op, cr)] = lab; print(r)
for band, k, op, cr in [("n_w50_toa", 2, 3, 0), ("n_w50_toa", 2, 3, 10), ("n_w30_toa", 2, 3, 5), ("n_w50_sr", 2, 3, 5)]:
    r, lab = run(band, k, op, cr); rows.append(r); labs[(band, k, op, cr)] = lab; print(r)
df = pd.DataFrame(rows); df.to_csv(os.path.join(OUT, f"08_outline_knobs_{TILE}_{YEAR}.csv"), index=False)
with open(os.path.join(OUT, f"08_outline_knobs_{TILE}_{YEAR}.txt"), "w") as f:
    f.write(f"Outline knobs, tile {TILE} {YEAR}, 10 m, every S2 scene, no scene cloud filter — {time.strftime('%Y-%m-%d %H:%M %Z')}\n")
    f.write(f"n_obs per pixel: median {np.median(n_obs):.0f}, 10/90 % {np.percentile(n_obs,10):.0f}/{np.percentile(n_obs,90):.0f}; TOA scenes {meta['n_toa_scenes']}, SR {meta['n_sr_scenes']}\n")
    f.write(f"Dunmire lakes inside tile: {len(dun)} ({dun_px.sum()*1e-4:.1f} km2). matched = >=50 % of the lake covered by a kept component; extra = kept components touching no Dunmire lake; chain = Dunmire lakes under one component\n\n")
    f.write(df.to_string(index=False) + "\n")

# figure: three views for k=1 vs k=2 (opening on), Dunmire in orange
views = {"crescent CW_0381": (-172_900, -2_144_700, 1_500), "dumbbell CW2018_1270": (-185_200, -2_109_500, 2_200),
         "wide view, tile centre": (meta["x0"] + 50_000, meta["y0"] + 50_000, 12_000)}
fig, axes = plt.subplots(len(views), 3, figsize=(17, 5.4 * len(views)))
for i, (name, (cx, cy, half)) in enumerate(views.items()):
    c0, r0 = ~tr * (cx - half, cy + half); c1, r1 = ~tr * (cx + half, cy - half)
    r0, r1, c0, c1 = int(max(r0, 0)), int(min(r1, N)), int(max(c0, 0)), int(min(c1, N)); ext = (cx - half, cx + half, cy - half, cy + half)
    sub = dun[dun.geometry.intersects(box(cx - half, cy - half, cx + half, cy + half))]
    for j, (ttl, arr, cmap, vmax) in enumerate([("scenes with NDWI>0.5 (count)", bands["n_w50_toa"][r0:r1, c0:c1], "viridis", None),
                                                ("k=1, opening 3 px, closing 5 px", (labs[("n_w50_toa", 1, 3, 5)][r0:r1, c0:c1] > 0), "Blues", 1),
                                                ("k=2, opening 3 px, closing 5 px", (labs[("n_w50_toa", 2, 3, 5)][r0:r1, c0:c1] > 0), "Blues", 1)]):
        ax = axes[i, j]; im = ax.imshow(arr, extent=ext, cmap=cmap, vmin=0, vmax=vmax if vmax else np.percentile(arr, 99) + 1, interpolation="nearest")
        if j == 0: plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        for g in sub.geometry:
            xx, yy = g.exterior.xy; ax.plot(xx, yy, color="tab:orange", lw=1.0)
        ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_title(f"{name}: {ttl}", fontsize=9); ax.set_xticks([]); ax.set_yticks([])
fig.suptitle(f"Sentinel-2 {YEAR}, tile {TILE}, 10 m, every scene, no scene cloud filter; orange = Dunmire {YEAR} outlines", fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(OUT, f"08_outline_knobs_{TILE}_{YEAR}.png"), dpi=100)
print("wrote", f"out/08_outline_knobs_{TILE}_{YEAR}.{{txt,csv,png}}")
