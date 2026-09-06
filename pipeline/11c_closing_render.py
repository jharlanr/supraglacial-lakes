"""Render the closing-accuracy chip strips from 11b's raw windows: panel 0 = ArcticDEM hillshade with the site polygon(s)
in red and per-year outlines in colour; panels 1-4 = true colour on the chosen scenes with the site polygon(s) dashed red.
Writes out/11_chips/{case_id}.png and out/11_chips/manifest.json.   Run:  TILE=19_39 nice -n 15 $(cat .python_env) scripts/11c_closing_render.py"""
import os, json, glob, numpy as np, scipy.ndimage as ndi, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LightSource; from affine import Affine; from shapely import wkt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(ROOT, "out"); CH = os.path.join(OUT, "11_chips"); RAW = os.path.join(CH, "raw")
TILE = os.environ.get("TILE", "19_39"); FORCE = os.environ.get("FORCE") == "1"
pref = f"08_s2counts_{TILE}_"; years = sorted(int(os.path.basename(f)[len(pref):].split("_")[0]) for f in glob.glob(os.path.join(OUT, pref + "*_meta.json")))
meta = json.load(open(os.path.join(OUT, f"{pref}{years[0]}_meta.json"))); tr = Affine(*meta["transform"]); N = meta["shape"][0]
cnts = {y: np.load(os.path.join(OUT, f"{pref}{y}_n_w50_toa.npy"), mmap_mode="r") for y in years}
dem = np.load(os.path.join(OUT, f"{TILE}_dem.npy")); dtr = Affine(*json.load(open(os.path.join(OUT, f"{TILE}_dem_meta.json")))["transform"])
polys = {k: wkt.loads(v) for k, v in json.load(open(os.path.join(OUT, f"11_sites3413_{TILE}.json"))).items()}
def disk(r):
    y, x = np.ogrid[-r:r + 1, -r:r + 1]; return (x * x + y * y) <= r * r
def draw_polys(ax, sids, **kw):
    for s in sids:
        g = polys[s]; parts = g.geoms if hasattr(g, "geoms") else [g]
        for p in parts:
            ax.plot(*p.exterior.xy, **kw)
            for h in p.interiors: ax.plot(*h.xy, **kw)
ls = LightSource(315, 45); cmap = plt.get_cmap("viridis", len(years)); mp = os.path.join(CH, "manifest.json"); manifest = json.load(open(mp)) if os.path.exists(mp) else {}
for jp in sorted(glob.glob(os.path.join(RAW, "*.json"))):
    c = json.load(open(jp)); cid = c["case_id"]; png = os.path.join(CH, f"{cid}.png")
    if not FORCE and os.path.exists(png) and cid in manifest and manifest[cid].get("tile") == TILE: continue
    if any(str(d).startswith("ERROR") for d in c["dates"].values()) and not all(os.path.exists(os.path.join(RAW, f"{cid}_{y}.npy")) for y in c["seasons"]): pass
    wx0, wy0, wx1, wy1 = c["window"]; ext = (wx0, wx1, wy0, wy1); ns = len(c["seasons"])
    fig, axes = plt.subplots(1, 1 + ns, figsize=(4.2 * (1 + ns), 4.7))
    c0, r0 = ~tr * (wx0, wy1); c1, r1 = ~tr * (wx1, wy0); r0, r1, c0, c1 = int(max(r0, 0)), int(min(r1, N)), int(max(c0, 0)), int(min(c1, N))
    dc0, dr0 = ~dtr * (wx0, wy1); dc1, dr1 = ~dtr * (wx1, wy0); dr0, dr1, dc0, dc1 = int(max(dr0, 0)), int(dr1) + 1, int(max(dc0, 0)), int(dc1) + 1
    dw = dem[dr0:dr1, dc0:dc1]; dext = ((dtr * (dc0, 0))[0], (dtr * (dc1, 0))[0], (dtr * (0, dr1))[1], (dtr * (0, dr0))[1])
    axes[0].imshow(ls.hillshade(np.nan_to_num(dw, nan=np.nanmean(dw)), vert_exag=3, dx=32, dy=32), cmap="gray", extent=dext)
    cext = ((tr * (c0, 0))[0], (tr * (c1, 0))[0], (tr * (0, r1))[1], (tr * (0, r0))[1])
    for k, y in enumerate(years):
        wmask = ndi.binary_closing(np.asarray(cnts[y][r0:r1, c0:c1]) >= 1, structure=disk(5))
        if wmask.any(): axes[0].contour(np.flipud(wmask), levels=[0.5], colors=[cmap(k)], linewidths=0.7, extent=cext, origin="lower")
    draw_polys(axes[0], c["site_ids"], color="red", linewidth=1.6)
    head = f"{cid}   [{c['pool']}, {c['n_pieces']} pieces" + (f", gap {c['gap_m']} m" if c["gap_m"] is not None else "") + "]"
    axes[0].set_title(f"DEM hillshade; red = site polygon(s); colours = outlines {years[0]}-{years[-1]}", fontsize=8)
    for ax, y in zip(axes[1:], c["seasons"]):
        f = os.path.join(RAW, f"{cid}_{y}.npy")
        if os.path.exists(f): ax.imshow(np.load(f), extent=ext); ax.set_title(f"{c['dates'][str(y)]}  Sentinel-2 true colour", fontsize=9)
        else: ax.text(0.5, 0.5, f"{y}: {c['dates'].get(str(y), 'missing')}", ha="center", va="center", transform=ax.transAxes, fontsize=7)
        draw_polys(ax, c["site_ids"], color="red", linewidth=0.9, linestyle="--")
    for ax in axes: ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_xticks([]); ax.set_yticks([])
    axes[1].plot([wx0 + 80, wx0 + 580], [wy0 + 80, wy0 + 80], color="white", linewidth=3); axes[1].text(wx0 + 330, wy0 + 130, "500 m", color="white", ha="center", fontsize=8)
    fig.suptitle(head, fontsize=10, y=0.995); fig.tight_layout(); fig.savefig(png, dpi=85); plt.close(fig)
    manifest[cid] = dict(tile=TILE, pool=c["pool"], site_ids=c["site_ids"], n_pieces=c["n_pieces"], gap_m=c["gap_m"], seasons=c["seasons"], dates=c["dates"], side_m=wx1 - wx0, png=os.path.basename(png))
json.dump(manifest, open(mp, "w"), indent=1); print("rendered", len(manifest), "chips")
