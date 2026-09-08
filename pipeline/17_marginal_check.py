"""Eyeball check for the ice-edge decision (spec 13e.1): the eleven sites on 29_45 that the OLD
"any non-ice contact" rule deleted, shown on the persistence surface with the ice-sheet exterior and
the nunataks drawn in.  Under the new exterior rule four are still excluded and seven come back --- this
figure is so Josh can say whether that is the right cut.

Panels are titled EXCLUDED / KEPT, with area, the seasons the site held water, and whether it touches
the ice-sheet exterior or only a nunatak.

Run:  nice -n 15 $(cat .python_env) scripts/17_marginal_check.py     (TILE=29_45)
Writes out/17_marginal_check_{TILE}.png and out/17_marginal_check_{TILE}.txt
"""
import os, json, numpy as np, pandas as pd, geopandas as gpd, scipy.ndimage as ndi
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from rasterio.features import rasterize
from affine import Affine
from pyproj import Transformer

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(ROOT, "out")
TILE = os.environ.get("TILE", "29_45")

# the eleven the old rule deleted, as (serial, lat, lon, km2) read off their pre-rebuild IDs
DELETED = [("G00099", 78.8512, -20.6892, 9.65), ("G00147", 78.7931, -20.8414, 3.16),
           ("G00153", 78.7857, -20.9957, 2.97), ("G00218", 78.6508, -23.0962, 1.47),
           ("G00244", 78.5486, -22.8712, 1.11), ("G00257", 78.4894, -22.3648, 0.68),
           ("G00273", 78.4025, -23.2577, 0.44), ("G00247", 78.5397, -23.2952, 0.36),
           ("G00230", 78.6013, -23.0057, 0.35), ("G00315", 78.2016, -23.6620, 0.25),
           ("G00318", 78.1697, -23.2427, 0.09)]

meta = json.load(open(os.path.join(OUT, f"08_s2counts_{TILE}_2019_meta.json")))
tr = Affine(*meta["transform"]); N = meta["shape"][0]; x0, y1 = meta["x0"], meta["y1"]
nsea = np.load(os.path.join(OUT, f"16_nsea_{TILE}.npy"))
years = 10

im = np.load(os.path.join(OUT, f"{TILE}_icemask_150m.npy"))
ex = np.load(os.path.join(OUT, f"{TILE}_extnonice_150m.npy"))
imm = json.load(open(os.path.join(OUT, f"{TILE}_icemask_150m_meta.json")))
xc = x0 + 10 * (np.arange(N) + 0.5); yc = y1 - 10 * (np.arange(N) + 0.5)
rows = np.clip(np.rint((imm["y_row0"] - yc) / 150).astype(int), 0, im.shape[0] - 1)
cols = np.clip(np.rint((xc - imm["x_col0"]) / 150).astype(int), 0, im.shape[1] - 1)
ice = im[rows[:, None], cols[None, :]]; extn = ex[rows[:, None], cols[None, :]]

exc = pd.read_csv(os.path.join(OUT, f"09_excluded_{TILE}.csv"))   # written by 08e: what the rule dropped
EXC_XY = list(zip(exc.x3413, exc.y3413))
sites = gpd.read_file(os.path.join(OUT, f"09_sites_{TILE}.geojson")).to_crs("EPSG:3413")
sr = rasterize(((g, 1) for g in sites.geometry), out_shape=(N, N), transform=tr, dtype=np.uint8)
yrs = pd.read_csv(os.path.join(OUT, f"09_site_years_{TILE}.csv"))
acols = [c for c in yrs.columns if c.startswith("a2")]
yrs[acols] = yrs[acols].apply(pd.to_numeric, errors="coerce").fillna(0)

tf = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
lines = [f"Ice-edge decision check, tile {TILE} — the 11 sites the old 'any non-ice contact' rule deleted",
         "KEPT/EXCLUDED as 08e actually decided (out/09_excluded_{TILE}.csv), not inferred", ""]

fig, axes = plt.subplots(3, 4, figsize=(15.5, 12), squeeze=False)
for ax, (name, lat, lon, km2) in zip(axes.ravel(), DELETED):
    X, Y = tf.transform(lon, lat)
    c = int((X - x0) / 10); r = int((y1 - Y) / 10)
    pad = max(120, int(np.sqrt(km2 * 1e6) / 10 * 1.5))
    r0, r1 = max(0, r - pad), min(N, r + pad); c0, c1 = max(0, c - pad), min(N, c + pad)
    v = nsea[r0:r1, c0:c1].astype(float)
    kept = sr[r0:r1, c0:c1].astype(bool)
    inside = not any(abs(X - ex_) < 400 and abs(Y - ey_) < 400 for ex_, ey_ in EXC_XY)   # excluded if 08e dropped it
    loc = nsea[max(0, r - 40):r + 40, max(0, c - 40):c + 40]
    peak = int(loc.max()) if loc.size else 0
    # background: nunatak grey, exterior non-ice dark
    bg = np.zeros(v.shape); bg[~ice[r0:r1, c0:c1] & ~extn[r0:r1, c0:c1]] = 1; bg[extn[r0:r1, c0:c1]] = 2
    ax.imshow(bg, cmap=ListedColormap(["#f2f2f2", "#b8a888", "#4a5a7a"]), vmin=0, vmax=2, interpolation="nearest")
    cm = plt.get_cmap("viridis", years + 1)
    ax.imshow(np.ma.masked_where(v == 0, v), cmap=cm, norm=BoundaryNorm(np.arange(-0.5, years + 1.5), cm.N), interpolation="nearest")
    if kept.any(): ax.contour(kept.astype(float), levels=[0.5], colors="r", linewidths=0.9)
    ax.set_title(f"{name}  {km2:.2f} km²  —  {'KEPT' if inside else 'EXCLUDED'}\n"
                 f"wettest pixel held water in {peak} of 10 seasons", fontsize=9,
                 color=("#14532d" if inside else "#7f1d1d"))
    ax.set_xticks([]); ax.set_yticks([])
    lines.append(f"{name}  {km2:6.2f} km2  {'KEPT    ' if inside else 'EXCLUDED'}  "
                 f"wettest pixel water in {peak}/10 seasons")
axes.ravel()[-1].axis("off")
axes.ravel()[-1].text(0.02, 0.6,
    "persistence colour = seasons the pixel held water (1–10)\n"
    "red outline = the site as the registry now keeps it\n"
    "grey = nunatak (non-ice enclosed by ice)\n"
    "blue-grey = ice-sheet exterior (ocean / coastal land)\n\n"
    "The old rule deleted all eleven because the 50 m closing\n"
    "bridged a non-ice pixel. The new rule asks WHICH non-ice.",
    fontsize=9, va="top")
f = os.path.join(OUT, f"17_marginal_check_{TILE}.png")
fig.suptitle(f"Tile {TILE}: what the ice-edge rule change does (spec §13e.1)", fontsize=12)
fig.tight_layout(); fig.savefig(f, dpi=125, bbox_inches="tight"); plt.close(fig)
open(os.path.join(OUT, f"17_marginal_check_{TILE}.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines)); print(f"\nfigure: {f}")
