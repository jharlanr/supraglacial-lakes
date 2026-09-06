"""Cases for the closing-accuracy test (spec section 3): every place where the 150 m site closing had to decide.
Pool A ("joined"): sites made of >= 2 separate pieces of the ten-season union before the closing.
Pool B ("apart"):  pairs of distinct sites whose polygons are within NEAR_M (default 1000 m) of each other.
Pool C ("lid"):    sites (not in A) that held >= 2 separate outlines in the same season at least once (ice lid, or two lakes?).
Writes out/11_cases_{TILE}.csv (case_id, pool, site_ids, n_pieces, gap_m, window x/y bounds in EPSG:3413, best seasons).
Run:  TILE=19_39 nice -n 15 $(cat .python_env) scripts/11a_closing_cases.py"""
import os, json, glob, numpy as np, pandas as pd, geopandas as gpd, scipy.ndimage as ndi
from rasterio.features import rasterize; from affine import Affine; from shapely.geometry import box
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(ROOT, "out")
TILE = os.environ.get("TILE", "19_39"); NEAR_M = float(os.environ.get("NEAR_M", "1000")); MIN_PX = 500; FILL = 0.5; S8 = np.ones((3, 3), bool)
def disk(r):
    y, x = np.ogrid[-r:r + 1, -r:r + 1]; return (x * x + y * y) <= r * r
pref = f"08_s2counts_{TILE}_"; years = sorted(int(os.path.basename(f)[len(pref):].split("_")[0]) for f in glob.glob(os.path.join(OUT, pref + "*_meta.json")))
meta = json.load(open(os.path.join(OUT, f"{pref}{years[0]}_meta.json"))); tr = Affine(*meta["transform"]); N = meta["shape"][0]; x0, y0, x1, y1 = meta["x0"], meta["y0"], meta["x1"], meta["y1"]
im = np.load(os.path.join(OUT, f"{TILE}_icemask_150m.npy")); imm = json.load(open(os.path.join(OUT, f"{TILE}_icemask_150m_meta.json")))
xc = x0 + 10 * (np.arange(N) + 0.5); yc = y1 - 10 * (np.arange(N) + 0.5)
ice = im[np.clip(np.rint((imm["y_row0"] - yc) / 150).astype(int), 0, im.shape[0] - 1)[:, None], np.clip(np.rint((xc - imm["x_col0"]) / 150).astype(int), 0, im.shape[1] - 1)[None, :]]
union = np.zeros((N, N), bool); ylabs = {}
for y in years:  # same per-year outline recipe as 08e
    arr = np.where(ice, np.load(os.path.join(OUT, f"{pref}{y}_n_w50_toa.npy")), 0)
    w = ndi.binary_closing(arr >= 1, structure=disk(5)); lab, _ = ndi.label(w, structure=S8)
    cnt = np.bincount(lab.ravel()); ids = np.flatnonzero(cnt >= MIN_PX); ids = ids[ids > 0]
    keep = ids[ndi.mean((arr >= 1).astype(float), lab, ids) >= FILL]
    keep = keep[ndi.maximum(ndi.binary_erosion(np.isin(lab, keep), structure=disk(2)).astype(np.uint8), lab, keep).astype(bool)]  # core test, as in 08e
    union |= np.isin(lab, keep); ylabs[y] = np.where(np.isin(lab, keep), lab, 0)
pieces, npc = ndi.label(union, structure=S8); print(f"union pieces before closing: {npc}")
multi = {}  # per site: max number of separate outlines in one season
gdf = gpd.read_file(os.path.join(OUT, f"09_sites_{TILE}.geojson")).to_crs(3413)
sy = pd.read_csv(os.path.join(OUT, f"09_site_years_{TILE}.csv"))
site_r = rasterize(((g, int(s)) for g, s in zip(gdf.geometry, gdf.serial + 1)), out_shape=(N, N), transform=tr, fill=0, dtype="int32")
# pieces per site: a piece belongs to the site covering most of it
pc = ndi.maximum(site_r, pieces, index=np.arange(1, npc + 1)); pieces_per_site = pd.Series(pc[pc > 0]).value_counts()
def best_seasons(sid, k=4):
    s = sy[sy.site_id == sid] if "site_id" in sy else sy[sy.site == sid]
    cols = [c for c in s.columns if c.startswith("area_") or c.isdigit()]
    return ""
print("site_years columns:", sy.columns.tolist()[:8])
rows = []
for serial1, n in pieces_per_site.items():
    if n < 2: continue
    r = gdf[gdf.serial == serial1 - 1].iloc[0]; b = r.geometry.bounds
    rows.append(dict(case_id=f"A_{r.site_id}", pool="joined", site_ids=r.site_id, n_pieces=int(n), gap_m=np.nan, xmin=b[0], ymin=b[1], xmax=b[2], ymax=b[3], area_km2=r.area_km2, n_years=r.n_years))
inA = set(c["site_ids"] for c in rows)
for y in years:
    yl = ylabs[y]; m = yl > 0; owner = ndi.maximum(site_r, yl, index=np.arange(1, yl.max() + 1)); vc = pd.Series(owner[owner > 0]).value_counts()
    for serial1, n in vc.items():
        if n >= 2: multi[serial1] = max(multi.get(serial1, 0), int(n))
for serial1, n in multi.items():
    r = gdf[gdf.serial == serial1 - 1].iloc[0]
    if r.site_id in inA: continue
    b = r.geometry.bounds; rows.append(dict(case_id=f"C_{r.site_id}", pool="lid", site_ids=r.site_id, n_pieces=n, gap_m=np.nan, xmin=b[0], ymin=b[1], xmax=b[2], ymax=b[3], area_km2=r.area_km2, n_years=r.n_years))
sidx = gdf.sindex
for i, r in gdf.iterrows():
    for j in sidx.query(r.geometry.buffer(NEAR_M), predicate="intersects"):
        if j <= i: continue
        q = gdf.iloc[j]; d = r.geometry.distance(q.geometry)
        if 0 < d <= NEAR_M:
            b = r.geometry.union(q.geometry).bounds
            rows.append(dict(case_id=f"B_{r.site_id}__{q.site_id}", pool="apart", site_ids=f"{r.site_id};{q.site_id}", n_pieces=2, gap_m=round(d), xmin=b[0], ymin=b[1], xmax=b[2], ymax=b[3], area_km2=round(r.area_km2 + q.area_km2, 3), n_years=max(r.n_years, q.n_years)))
c = pd.DataFrame(rows); c.to_csv(os.path.join(OUT, f"11_cases_{TILE}.csv"), index=False)
print(c.pool.value_counts().to_dict()); print(c.groupby("pool").n_pieces.describe()[["mean", "max"]]); print("gap<=400:", int((c.gap_m <= 400).sum())); print("gap m quantiles", c[c.pool == "apart"].gap_m.quantile([.25, .5, .75]).round().tolist())
print("window size m: median", round(((c.xmax - c.xmin).clip(lower=c.ymax - c.ymin)).median()), "max", round(((c.xmax - c.xmin).clip(lower=c.ymax - c.ymin)).max()))
