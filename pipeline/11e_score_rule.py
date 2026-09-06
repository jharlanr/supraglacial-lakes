"""Score the appendage rule (spec 8b) against both label sets on the 105 closing-test cases, for several parameter settings,
without rebuilding the registry: rebuilds the ten-season union pieces (as 08e), applies the rule per setting, and asks for each
labelled case whether its bodies end up in one site. Cases labelled 3 (can't tell) are excluded.
Writes out/11_rule_scores_{TILE}.txt.   Run:  TILE=19_39 nice -n 15 $(cat .python_env) scripts/11e_score_rule.py"""
import os, json, glob, numpy as np, pandas as pd, geopandas as gpd, scipy.ndimage as ndi
from rasterio.features import shapes; from affine import Affine; from shapely.geometry import shape; from shapely.ops import unary_union; from shapely import wkt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(ROOT, "out"); TILE = os.environ.get("TILE", "19_39")
SETTINGS = [(50, 150, 0.2), (50, 200, 0.2), (50, 200, 0.25), (50, 250, 0.2), (150, 150, 1.0)]  # (JOIN_ALL_M, JOIN_APP_M, APP_RATIO); last = "join everything within 150 m"
MIN_PX = 500; FILL = 0.5; S8 = np.ones((3, 3), bool)
def disk(r):
    y, x = np.ogrid[-r:r + 1, -r:r + 1]; return (x * x + y * y) <= r * r
pref = f"08_s2counts_{TILE}_"; years = sorted(int(os.path.basename(f)[len(pref):].split("_")[0]) for f in glob.glob(os.path.join(OUT, pref + "*_meta.json")))
meta = json.load(open(os.path.join(OUT, f"{pref}{years[0]}_meta.json"))); tr = Affine(*meta["transform"]); N = meta["shape"][0]; x0, y0, x1, y1 = meta["x0"], meta["y0"], meta["x1"], meta["y1"]
im = np.load(os.path.join(OUT, f"{TILE}_icemask_150m.npy")); imm = json.load(open(os.path.join(OUT, f"{TILE}_icemask_150m_meta.json")))
xc = x0 + 10 * (np.arange(N) + 0.5); yc = y1 - 10 * (np.arange(N) + 0.5)
ice = im[np.clip(np.rint((imm["y_row0"] - yc) / 150).astype(int), 0, im.shape[0] - 1)[:, None], np.clip(np.rint((xc - imm["x_col0"]) / 150).astype(int), 0, im.shape[1] - 1)[None, :]]
union = np.zeros((N, N), bool)
for y in years:
    arr = np.where(ice, np.load(os.path.join(OUT, f"{pref}{y}_n_w50_toa.npy")), 0); w = ndi.binary_closing(arr >= 1, structure=disk(5)); lab, _ = ndi.label(w, structure=S8)
    cnt = np.bincount(lab.ravel()); ids = np.flatnonzero(cnt >= MIN_PX); ids = ids[ids > 0]; keep = ids[ndi.mean((arr >= 1).astype(float), lab, ids) >= FILL]
    keep = keep[ndi.maximum(ndi.binary_erosion(np.isin(lab, keep), structure=disk(2)).astype(np.uint8), lab, keep).astype(bool)]; union |= np.isin(lab, keep)
plab, npieces = ndi.label(union, structure=S8); ppx = np.bincount(plab.ravel())
pgeom = {}
for geom, val in shapes(plab.astype(np.int32), mask=plab > 0, transform=tr): pgeom.setdefault(int(val), []).append(shape(geom))
pg = gpd.GeoDataFrame({"piece": list(pgeom)}, geometry=[unary_union(v) for v in pgeom.values()], crs=3413); pg["px"] = ppx[pg.piece.values]; sidx = pg.sindex
# candidate pairs once (within 300 m)
pairs = []
for i, r in pg.iterrows():
    for j in sidx.query(r.geometry.buffer(300), predicate="intersects"):
        q = pg.iloc[j]
        if int(q.piece) <= int(r.piece): continue
        pairs.append((int(r.piece), int(q.piece), r.geometry.distance(q.geometry), min(r.px, q.px) / max(r.px, q.px)))
def roots(all_m, app_m, ratio_max):
    parent = {int(p): int(p) for p in pg.piece}
    def find(a):
        while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for a, b, d, ratio in pairs:
        if d <= all_m or (d <= app_m and ratio < ratio_max):
            ra, rb = find(a), find(b)
            if ra != rb: parent[rb] = ra
    return {p: find(p) for p in parent}
# labelled cases -> the pieces each old site polygon covers (old = the 297-site build the chips were drawn from)
polys = json.load(open(os.path.join(OUT, f"11_sites3413_{TILE}.json"))); both = pd.read_csv(os.path.join(OUT, f"11_labels_both_{TILE}.csv"))
def pieces_of(sid):
    g = wkt.loads(polys[sid]); return [int(pg.iloc[j].piece) for j in sidx.query(g, predicate="intersects") if pg.iloc[j].geometry.intersection(g).area > 0.5 * pg.iloc[j].geometry.area]
case_pieces = {}
for _, r in both.iterrows():
    sids = r.case_id[2:].split("__"); case_pieces[r.case_id] = [pieces_of(s) for s in sids]
lines = [f"Appendage-rule scoring on tile {TILE}, {len(both)} labelled cases (label 3 excluded per labeller). 'joined' = the case's bodies end in one site."]
for all_m, app_m, rmax in SETTINGS:
    rt = roots(all_m, app_m, rmax); rows = []
    for _, r in both.iterrows():
        groups = case_pieces[r.case_id]; allp = [p for g in groups for p in g]
        if r.pool == "apart":
            joined = len(groups) == 2 and len(groups[0]) and len(groups[1]) and any(rt[a] == rt[b] for a in groups[0] for b in groups[1])
        else:
            joined = len(set(rt[p] for p in allp)) <= 1 if allp else True
        rows.append(dict(case_id=r.case_id, pool=r.pool, joined=bool(joined), josh=r.josh, claude=r.claude))
    d = pd.DataFrame(rows); s = f"\n[{all_m} m always | {app_m} m if ratio < {rmax}]  sites joined: joined-pool {int(d[d.pool=='joined'].joined.sum())}/16, lid-pool {int(d[d.pool=='lid'].joined.sum())}/35, apart-pool {int(d[d.pool=='apart'].joined.sum())}/54"
    for who in ("josh", "claude"):
        v = d[d[who] != 3]; want = v[who] == 1; ok = (v.joined == want); fm = int(((~want) & v.joined).sum()); fs = int((want & ~v.joined).sum())
        s += f"\n   vs {who}: agree {int(ok.sum())}/{len(v)} ({ok.mean()*100:.0f} %); false merges {fm}, false splits {fs}"
        if fm + fs: s += "  -> " + "; ".join(f"{c.case_id[2:8]}{'..' if c.pool=='apart' else ''}:{'FM' if c.joined else 'FS'}" for c in v[~ok].itertuples())
    lines.append(s); print(s)
open(os.path.join(OUT, f"11_rule_scores_{TILE}.txt"), "w").write("\n".join(lines) + "\n")
