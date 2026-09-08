"""Ten-season lake sites + registry for one tile (v0.0 recipe, notes/spec_v0.0_lake_sites.md):
per-year outlines from the 08b counts (TOA NDWI>0.5, any scene, closing 5 px, >=0.05 km2, fill>=0.5) ->
union over seasons -> pieces joined by the appendage rule (gap <= 50 m always; gap <= 250 m if the smaller is < 1/5 of the larger; spec 8b) -> sites ->
registry (ID = centroid lat/lon at registration), per-year presence and area, DEM attributes (32 m sinks /
sub-basins / depth from 04b, reported not enforced), crosswalk to Dunmire 2018/2019, figures.
Outputs: out/09_sites_{TILE}.geojson, out/09_registry_{TILE}.csv, out/09_site_years_{TILE}.csv,
         out/09_dunmire_crosswalk_{TILE}.csv, out/09_sites_{TILE}.txt, out/09_sites_{TILE}_{map,showcase,stats}.png
Run:  nice -n 15 $(cat .python_env) scripts/08e_sites_registry.py        (TILE=19_39; JOIN_ALL_M=50 JOIN_APP_M=250 APP_RATIO=0.2 MINW_PX=0 EXCLUDE_TOUCHING=1)
"""
import os, json, glob, time
import numpy as np, pandas as pd, geopandas as gpd
from scipy import ndimage as ndi
from rasterio.features import rasterize, shapes
from rasterio.transform import Affine
from shapely.geometry import shape, box
from shapely.ops import unary_union
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(ROOT, "out")
TILE = os.environ.get("TILE", "19_39"); JOIN_ALL_M = float(os.environ.get("JOIN_ALL_M", "50")); JOIN_APP_M = float(os.environ.get("JOIN_APP_M", "250")); APP_RATIO = float(os.environ.get("APP_RATIO", "0.2")); FILL = 0.5; MINW_PX = int(os.environ.get("MINW_PX", "0")); EXCLUDE_TOUCHING = os.environ.get("EXCLUDE_TOUCHING", "0") == "1"   # Josh 2026-09-07: no ice-marginal exclusion; on ice is enough
# spec 13a: the definition is stated in METRES and every length is rounded UP to whole pixels, so the
# same rule gives the same lake at 10 m and at 20 m.  TOUCH_RULE: "core" (spec 13c.1, the site's 50 m
# core must reach non-ice) or "any" (the old rule: one fringe pixel of contact excluded the site).
CLOSE_M = float(os.environ.get("CLOSE_M", "50")); CORE_M = float(os.environ.get("CORE_M", "50"))
MIN_AREA_KM2 = float(os.environ.get("MIN_AREA_KM2", "0.05")); TOUCH_RULE = os.environ.get("TOUCH_RULE", "exterior")
pref = f"08_s2counts_{TILE}_"
years = sorted(int(os.path.basename(f)[len(pref):].split("_")[0]) for f in glob.glob(os.path.join(OUT, pref + "*_meta.json")))
meta = json.load(open(os.path.join(OUT, f"{pref}{years[0]}_meta.json"))); tr = Affine(*meta["transform"]); N = meta["shape"][0]
PX = abs(tr.a)                                    # pixel size in m, from the export (10 or 20)
import math
R_CLOSE = max(1, math.ceil(CLOSE_M / PX))          # closing disk radius, px  (50 m -> 5 px at 10 m, 3 px at 20 m)
R_CORE  = max(1, math.ceil(CORE_M / 2 / PX))       # erosion radius for the core, px (50 m across)
MIN_PX  = int(math.ceil(MIN_AREA_KM2 * 1e6 / PX ** 2))   # 0.05 km2 -> 500 px at 10 m, 125 px at 20 m
x0, y0, x1, y1 = meta["x0"], meta["y0"], meta["x1"], meta["y1"]
S8 = np.ones((3, 3), bool)
def disk(r):
    y, x = np.ogrid[-r:r + 1, -r:r + 1]; return (x * x + y * y) <= r * r
lines = [f"Sites v0.0, tile {TILE}, {PX:.0f} m, seasons {years}, closing {CLOSE_M:.0f} m ({R_CLOSE} px), core {CORE_M:.0f} m, min area {MIN_AREA_KM2} km2 ({MIN_PX} px), touch rule {TOUCH_RULE}, appendage {JOIN_ALL_M:.0f} m / {JOIN_APP_M:.0f} m / ratio {APP_RATIO} — {time.strftime('%Y-%m-%d %H:%M %Z')}"]
def say(s): print(s); lines.append(s)

# 0. domain: the ice sheet (BedMachine v6 mask, grounded + floating ice, from 08g), nearest-upsampled to the 10 m grid
imf = os.path.join(OUT, f"{TILE}_icemask_150m.npy")
if os.path.exists(imf):
    im = np.load(imf); imm = json.load(open(os.path.join(OUT, f"{TILE}_icemask_150m_meta.json")))
    xc = x0 + PX * (np.arange(N) + 0.5); yc = y1 - PX * (np.arange(N) + 0.5)
    cols = np.clip(np.rint((xc - imm["x_col0"]) / 150).astype(int), 0, im.shape[1] - 1); rows = np.clip(np.rint((imm["y_row0"] - yc) / 150).astype(int), 0, im.shape[0] - 1)
    ice = im[rows[:, None], cols[None, :]]; say(f"domain: BedMachine v6 ice mask, {ice.mean()*100:.1f} % of the tile is ice")
    exf = os.path.join(OUT, f"{TILE}_extnonice_150m.npy")   # non-ice connected to outside the ice sheet (08g)
    ext_nonice = np.load(exf)[rows[:, None], cols[None, :]] if os.path.exists(exf) else None
else:
    ice = np.ones((N, N), bool); ext_nonice = None; say("domain: no ice mask file (08g) — whole tile treated as ice")

# 1. per-year outlines
ylab = {}; yrows = []
for y in years:
    arr = np.load(os.path.join(OUT, f"{pref}{y}_n_w50_toa.npy")); arr = np.where(ice, arr, 0)
    # Josh 2026-09-07: a lake is water sitting ON ICE — being at the margin is fine.  The water mask is
    # already clipped to ice, so the only way a site ever held non-ice was rule 6's closing bridging across
    # it; re-applying the mask after the closing stops that, and then nothing needs excluding at all.
    w = ndi.binary_closing(arr >= 1, structure=disk(R_CLOSE)) & ice; lab, _ = ndi.label(w, structure=S8)
    cnt = np.bincount(lab.ravel()); ids = np.flatnonzero(cnt >= MIN_PX); ids = ids[ids > 0]
    fill = ndi.mean((arr >= 1).astype(float), lab, ids); keep = ids[fill >= FILL]
    core = ndi.maximum(ndi.binary_erosion(np.isin(lab, keep), structure=disk(R_CORE)).astype(np.uint8), lab, keep).astype(bool); keep = keep[core]  # must contain a 50 m wide core (drops swath-edge lines)
    M = np.isin(lab, keep)
    if MINW_PX:  # minimum width 2*MINW_PX*10 m (Josh 2026-09-06: 40 m): open, then grow back inside the outline so shapes keep, channels and necks go
        M = ndi.binary_dilation(ndi.binary_opening(M, structure=disk(MINW_PX)), structure=disk(MINW_PX)) & M
        lab, _ = ndi.label(M, structure=S8); cnt = np.bincount(lab.ravel()); keep = np.flatnonzero(cnt >= MIN_PX); keep = keep[keep > 0]
    ylab[y] = np.where(np.isin(lab, keep), lab, 0)
    m = json.load(open(os.path.join(OUT, f"{pref}{y}_meta.json")))
    yrows.append(dict(year=y, scenes=m["n_toa_scenes"], outlines=len(keep), area_km2=round(cnt[keep].sum() * 1e-4, 1)))
    say(f"  {y}: {m['n_toa_scenes']} scenes, {len(keep)} outlines, {cnt[keep].sum()*1e-4:.1f} km2")
pd.DataFrame(yrows).to_csv(os.path.join(OUT, f"09_year_outlines_{TILE}.csv"), index=False)

# 2. sites
union = np.zeros((N, N), bool)
for y in years: union |= ylab[y] > 0
# appendage rule (spec 8b, Josh 2026-09-05): union pieces are one site if their gap <= JOIN_ALL_M (ice-lid case), or if the gap
# is <= JOIN_APP_M and the smaller piece is under APP_RATIO of the larger (an appendage); otherwise separate sites. No site-level closing.
plab, npieces = ndi.label(union, structure=S8); ppx = np.bincount(plab.ravel())
pgeom = {}
for geom, val in shapes(plab.astype(np.int32), mask=plab > 0, transform=tr): pgeom.setdefault(int(val), []).append(shape(geom))
pg = gpd.GeoDataFrame({"piece": list(pgeom)}, geometry=[unary_union(v) for v in pgeom.values()], crs=3413); pg["px"] = ppx[pg.piece.values]
parent = {int(p): int(p) for p in pg.piece}
def find(a):
    while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
    return a
sidx = pg.sindex; n_all = n_app = 0; joins = []
for i, r in pg.iterrows():
    for j in sidx.query(r.geometry.buffer(JOIN_APP_M), predicate="intersects"):
        q = pg.iloc[j]
        if int(q.piece) <= int(r.piece): continue
        d = r.geometry.distance(q.geometry); ratio = min(r.px, q.px) / max(r.px, q.px)
        if d <= JOIN_ALL_M: kind = "lid"
        elif d <= JOIN_APP_M and ratio < APP_RATIO: kind = "appendage"
        else: continue
        a, b = find(int(r.piece)), find(int(q.piece))
        if a != b: parent[b] = a
        joins.append(dict(piece_a=int(r.piece), piece_b=int(q.piece), gap_m=round(d), ratio=round(ratio, 3), kind=kind)); n_all += kind == "lid"; n_app += kind == "appendage"
root = {int(p): find(int(p)) for p in pg.piece}; site_of_root = {rt: k + 1 for k, rt in enumerate(sorted(set(root.values())))}
lut = np.zeros(npieces + 1, np.int32)
for p, rt in root.items(): lut[p] = site_of_root[rt]
slab = lut[plab]; spx = np.bincount(slab.ravel()); ids = np.flatnonzero(spx >= MIN_PX); ids = ids[ids > 0]
slab = np.where(np.isin(slab, ids), slab, 0)
edge = ndi.distance_transform_edt(ice) * PX  # distance to non-ice (rock, ocean) on the BedMachine mask, m
if EXCLUDE_TOUCHING:  # water that reaches rock or ocean is an ice-marginal lake (How 2025's inventory), not a supraglacial one
    if TOUCH_RULE == "exterior":
        # spec 13c.1 (2026-09-07): a supraglacial lake may abut a NUNATAK; an ice-marginal lake reaches the
        # EDGE of the ice sheet.  The water mask is clipped to ice before the closing, so a site can only hold
        # non-ice pixels where rule 6's closing bridged them — under the old "any non-ice" rule a bridge over a
        # stray interior mask cell deleted the whole lake (three Dunmire lakes on 29_45).  Testing only against
        # non-ice connected to the outside of the ice sheet (08g, flood-filled from a padded window so it does
        # not depend on the tile cut) keeps those and still excludes the large marginal bodies.
        if ext_nonice is None: raise SystemExit("TOUCH_RULE=exterior needs out/{TILE}_extnonice_150m.npy — run 08g first")
        hit = ndi.maximum(ext_nonice.astype(np.uint8), slab, ids)
        e0 = np.where(np.atleast_1d(hit) > 0, 0.0, 1.0)      # 0 = reaches the ice-sheet exterior
    elif TOUCH_RULE == "core":
        # spec 13c.1: the site's 50 m CORE must reach non-ice.  The closing (rule 6) runs after the ice
        # mask, so it can push a fringe pixel across the boundary; under the old "any pixel" rule that
        # single pixel excluded the whole lake, which cost five Dunmire lakes on 29_45 sitting 335-541 m
        # inside the ice edge.  Using the core reuses rule 9's 50 m disk and adds no new constant.
        core = ndi.binary_erosion(slab > 0, structure=disk(R_CORE))
        e0 = ndi.minimum(np.where(core, edge, np.inf), slab, ids)
        e0 = np.where(np.isfinite(e0), e0, np.inf)      # a site with no core left: nothing to test, keep it
    else:
        e0 = ndi.minimum(edge, slab, ids)
    touching = ids[e0 == 0]
    if len(touching):   # record what was dropped, so the decision can be eyeballed (17_marginal_check.py)
        cy, cx = np.array(ndi.center_of_mass(slab > 0, slab, touching)).T
        pd.DataFrame(dict(excluded_label=touching, area_km2=(spx[touching] * PX ** 2 * 1e-6).round(4),
                          x3413=(x0 + (cx + 0.5) * PX).round(1), y3413=(y1 - (cy + 0.5) * PX).round(1),
                          rule=TOUCH_RULE)).to_csv(os.path.join(OUT, f"09_excluded_{TILE}.csv"), index=False)
    ids = ids[e0 > 0]; slab = np.where(np.isin(slab, ids), slab, 0)
    say(f"excluded {len(touching)} sites ({TOUCH_RULE} rule) reaching non-ice "
        f"({spx[touching].sum()*1e-4:.1f} km2); {len(ids)} sites remain")
pd.DataFrame(joins).to_csv(os.path.join(OUT, f"09_joins_{TILE}.csv"), index=False)
say(f"sites: {len(ids)} from {npieces} union pieces; {n_all} lid joins (gap <= {JOIN_ALL_M} m), {n_app} appendage joins (gap <= {JOIN_APP_M} m, ratio < {APP_RATIO}); total {spx[ids].sum()*1e-4:.1f} km2")

# 3. per-site per-year presence and water area (water = that year's outline pixels inside the site)
pres = {}
for y in years:
    a = np.bincount(slab[(slab > 0) & (ylab[y] > 0)], minlength=slab.max() + 1)
    pres[y] = a
site_years = pd.DataFrame({"site": ids, **{f"a{y}": np.round(pres[y][ids] * 1e-4, 4) for y in years}})
nyears = (site_years[[f"a{y}" for y in years]] > 0).sum(axis=1).values
first = [min((y for y in years if pres[y][i] > 0), default=-1) for i in ids]; last = [max((y for y in years if pres[y][i] > 0), default=-1) for i in ids]

# 4. DEM attributes (32 m products from 04a/04b), sampled at 10 m by nearest
dem_meta = json.load(open(os.path.join(OUT, f"{TILE}_dem_meta.json")))
def up(a32, smooth=False):  # 32 m tile array -> 10 m grid of this tile (same origin; both cover the 100 km tile)
    if smooth:  # bilinear, for the hillshade only
        z = ndi.zoom(np.nan_to_num(a32, nan=float(np.nanmin(a32))).astype("float32"), 3.2, order=1)
        out = np.empty((N, N), "float32"); h, w = min(N, z.shape[0]), min(N, z.shape[1]); out[:h, :w] = z[:h, :w]
        if h < N: out[h:, :] = out[h - 1:h, :]
        if w < N: out[:, w:] = out[:, w - 1:w]
        return out
    r = (np.arange(N) * PX // 32).astype(int).clip(0, a32.shape[0] - 1); c = (np.arange(N) * PX // 32).astype(int).clip(0, a32.shape[1] - 1)
    return a32[np.ix_(r, c)]
depth = up(np.load(os.path.join(OUT, f"{TILE}_depth.npy"))); sinks = up(np.load(os.path.join(OUT, f"{TILE}_sinks.npy"))); subs = up(np.load(os.path.join(OUT, f"{TILE}_subbasins.npy")))
frac_sink = ndi.mean((sinks > 0).astype(float), slab, ids); maxdepth = ndi.maximum(np.nan_to_num(depth), slab, ids)
main_sub = []
for i in ids:
    s = subs[slab == i]; s = s[s > 0]; main_sub.append(int(np.bincount(s).argmax()) if s.size else 0)
ice_edge_m = ndi.minimum(edge, slab, ids)
halfw = ndi.maximum(ndi.distance_transform_edt(slab > 0) * PX, slab, ids)  # widest point of the site (m); <= 60 m flags line-like sites for review  # distance from the site's nearest pixel to non-ice (BedMachine mask)
cy, cx = zip(*ndi.center_of_mass(slab > 0, slab, ids)); cxg = np.array([tr * (c, r) for r, c in zip(cy, cx)])
cent_in_sink = sinks[np.clip(np.array(cy).astype(int), 0, N - 1), np.clip(np.array(cx).astype(int), 0, N - 1)] > 0

# 5. polygons + IDs (centroid lat/lon, 4 decimals, EPSG:4326)
polys = {}
for geom, val in shapes(slab.astype("int32"), mask=slab > 0, transform=tr):
    polys.setdefault(int(val), []).append(shape(geom))
gdf = gpd.GeoDataFrame({"site": ids}, geometry=[unary_union(polys[i]) for i in ids], crs=3413)
elong = np.round((gdf.length ** 2 / (4 * np.pi * gdf.area)).values, 1)  # shape index: 1 = circle; > 15 flags line-like sites for review
ll = gpd.GeoSeries(gpd.points_from_xy(cxg[:, 0], cxg[:, 1]), crs=3413).to_crs(4326)
def mkid(serial, lat, lon):  # Josh's convention 2026-09-05: G<serial 5 digits>_<N|S><lat*1e4, 6 digits>_<E|W><lon*1e4, 7 digits>
    return f"G{serial:05d}_{'N' if lat >= 0 else 'S'}{round(abs(lat)*1e4):06d}_{'E' if lon >= 0 else 'W'}{round(abs(lon)*1e4):07d}"
# provisional serials for this tile: a north-to-south sweep (the Greenland-wide build numbers by drainage basin, then latitude)
order = np.argsort(-ll.y.values, kind="stable"); serial = np.empty(len(ids), int); serial[order] = np.arange(len(ids))
site_id = [mkid(int(n), p.y, p.x) for n, p in zip(serial, ll)]
seen = {}
for k, sid in enumerate(site_id):
    if sid in seen: seen[sid] += 1; site_id[k] = f"{sid}{seen[sid]}"
    else: seen[sid] = 1
reg = pd.DataFrame(dict(site_id=site_id, serial=serial, site=ids, lat=np.round(ll.y.values, 5), lon=np.round(ll.x.values, 5), x3413=np.round(cxg[:, 0]), y3413=np.round(cxg[:, 1]),
                        area_km2=np.round(spx[ids] * 1e-4, 4), n_years=nyears, first_seen=first, last_seen=last,
                        max_year_area_km2=np.round(np.max([pres[y][ids] for y in years], axis=0) * 1e-4, 4),
                        ice_edge_m=np.round(ice_edge_m), ice_marginal=ice_edge_m < 300, max_halfwidth_m=np.round(halfw), thin=halfw <= 60, elongation=elong, dem_frac_in_sink=np.round(frac_sink, 3), dem_centroid_in_sink=cent_in_sink, dem_max_depth_m=np.round(maxdepth, 1), dem_main_subbasin=main_sub,
                        registered="v0.0-" + time.strftime("%Y-%m-%d"), basis=f"S2 {years[0]}-{years[-1]}"))
gdf = gdf.merge(reg, on="site"); gdf.to_crs(4326).to_file(os.path.join(OUT, f"09_sites_{TILE}.geojson"), driver="GeoJSON")
reg.to_csv(os.path.join(OUT, f"09_registry_{TILE}.csv"), index=False); site_years.merge(reg[["site", "site_id"]]).to_csv(os.path.join(OUT, f"09_site_years_{TILE}.csv"), index=False)

# 6. Dunmire crosswalk
xw = []; tile_box = box(x0, y0, x1, y1)
for y in (2018, 2019):
    d = gpd.read_file(os.path.join(ROOT, "..", "labels", "dunmire", f"labels_{y}_volumes.geojson")).to_crs(3413); d = d[d.geometry.within(tile_box)]
    dl = rasterize(((g, i + 1) for i, g in enumerate(d.geometry)), out_shape=(N, N), transform=tr, dtype="int32"); dpx = np.bincount(dl.ravel())
    for i, r in enumerate(d.itertuples()):
        m = dl == i + 1; s = slab[m]; s = s[s > 0]
        if s.size == 0: xw.append(dict(dunmire_id=r.new_id, year=y, site_id="", frac=0.0, n_sites=0)); continue
        bc = np.bincount(s); top = bc.argmax(); n_sites = int((bc[1:] >= 0.1 * m.sum()).sum())
        xw.append(dict(dunmire_id=r.new_id, year=y, site_id=reg.set_index("site").loc[top, "site_id"], frac=round(bc[top] / m.sum(), 3), n_sites=n_sites))
xw = pd.DataFrame(xw); xw.to_csv(os.path.join(OUT, f"09_dunmire_crosswalk_{TILE}.csv"), index=False)
for y in (2018, 2019):
    g = xw[xw.year == y]; say(f"Dunmire {y}: {len(g)} lakes in tile; {int((g.frac>0).sum())} touch a site, {int((g.frac>=0.5).sum())} covered >=50 %, {int((g.n_sites>=2).sum())} span >=2 sites; sites holding >=2 {y} lakes: {int((g[g.frac>0].groupby('site_id').size()>=2).sum())}")

# 7. stats
say(f"ice-marginal sites (within 300 m of non-ice): {int((ice_edge_m < 300).sum())}; thin sites (widest point <= 60 m): {int((halfw <= 60).sum())}")
say(f"sites by years present (of {len(years)}): " + ", ".join(f"{k}:{v}" for k, v in sorted(pd.Series(nyears).value_counts().items())))
say(f"sites with centroid in a 32 m depression: {cent_in_sink.mean()*100:.0f} %; with >=50 % of area in one: {(frac_sink>=0.5).mean()*100:.0f} %; entirely outside: {(frac_sink==0).mean()*100:.0f} %")
for lo, hi in ((0, 0.1), (0.1, 0.2), (0.2, 0.5), (0.5, 1), (1, 100)):
    sel = (reg.area_km2 > lo) & (reg.area_km2 <= hi); say(f"  sites {lo}-{hi} km2: n={sel.sum()}, centroid in depression {cent_in_sink[sel].mean()*100:.0f} %, median years present {np.median(nyears[sel]) if sel.any() else 0:.0f}")
open(os.path.join(OUT, f"09_sites_{TILE}.txt"), "w").write("\n".join(lines) + "\n")

# 8. figures
dem = up(np.load(os.path.join(OUT, f"{TILE}_dem.npy")), smooth=True); ls = LightSource(315, 40)
fig, ax = plt.subplots(figsize=(13, 13)); hs = ls.hillshade(np.nan_to_num(dem, nan=np.nanmin(dem)), vert_exag=1, dx=10, dy=10)
ax.imshow(hs[::4, ::4], cmap="gray", extent=(x0, x1, y0, y1)); ax.imshow(np.where(sinks[::4, ::4] > 0, 1, np.nan), cmap="Blues_r", alpha=0.25, extent=(x0, x1, y0, y1), vmin=0, vmax=2)
gdf.plot(ax=ax, column="n_years", cmap="plasma", legend=True, legend_kwds={"label": "seasons with water", "shrink": 0.5}, edgecolor="k", linewidth=0.2)
ax.set_title(f"Tile {TILE}: {len(ids)} lake sites from the {years[0]}–{years[-1]} Sentinel-2 union; light blue = 32 m ArcticDEM depressions (attribute only)"); ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout(); fig.savefig(os.path.join(OUT, f"09_sites_{TILE}_map.png"), dpi=110); plt.close(fig)
views = {"crescent CW_0381": (-172_900, -2_144_700, 1_500), "dumbbell CW2018_1270": (-185_200, -2_109_500, 2_200)}
views = {k: v for k, v in views.items() if x0 <= v[0] <= x1 and y0 <= v[1] <= y1}
if len(views) < 2:  # other tiles: the two sites holding the most Dunmire lakes (ties by area), each framed to its own extent
    multi = xw[xw.site_id != ""].groupby("site_id").dunmire_id.nunique().sort_values(ascending=False)
    for sid in list(multi.index[:2]) + list(gdf.sort_values("area_km2", ascending=False).site_id):
        if len(views) >= 2: break
        g = gdf[gdf.site_id == sid].geometry.iloc[0]; b = g.bounds; views[f"{sid} ({multi.get(sid, 0)} Dunmire lakes)"] = ((b[0] + b[2]) / 2, (b[1] + b[3]) / 2, max(b[2] - b[0], b[3] - b[1]) * 0.8 + 400)
views["wide view"] = (x0 + 50_000, y0 + 50_000, 12_000)
cmap = plt.get_cmap("viridis", len(years)); fig, axes = plt.subplots(1, 3, figsize=(19, 6.5))
for ax, (vn, (cx_, cy_, half)) in zip(axes, views.items()):
    ext = (cx_ - half, cx_ + half, cy_ - half, cy_ + half); c0, r0 = ~tr * (ext[0], ext[3]); c1, r1 = ~tr * (ext[1], ext[2]); r0, r1, c0, c1 = int(max(r0, 0)), int(min(r1, N)), int(max(c0, 0)), int(min(c1, N))
    ax.imshow(ls.hillshade(np.nan_to_num(dem[r0:r1, c0:c1]), vert_exag=2, dx=10, dy=10), cmap="gray", extent=ext)
    sub = gdf[gdf.geometry.intersects(box(ext[0], ext[2], ext[1], ext[3]))]; sub.plot(ax=ax, facecolor="none", edgecolor="red", linewidth=1.6)
    for k, y in enumerate(years):
        yl = ylab[y][r0:r1, c0:c1]; ax.contour(np.flipud(yl > 0), levels=[0.5], colors=[cmap(k)], linewidths=0.8, extent=ext, origin="lower")
    for _, r in sub.iterrows(): ax.annotate(r.site_id if half < 5000 else r.site_id.split("_")[0], (r.geometry.centroid.x, r.geometry.centroid.y), fontsize=6 if half < 5000 else 5, ha="center", va="bottom", color="red", xytext=(0, 3), textcoords="offset points")
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(vn, fontsize=10)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(years[0] - 0.5, years[-1] + 0.5)); fig.colorbar(sm, ax=axes, shrink=0.6, label="season of the per-year outline")
fig.suptitle(f"Tile {TILE}: red = v0.0 site polygon (union of {years[0]}–{years[-1]} outlines + 150 m closing); coloured lines = per-year outlines; labels = Lake IDs (serial only in the wide view)", fontsize=10)
fig.savefig(os.path.join(OUT, f"09_sites_{TILE}_showcase.png"), dpi=100); plt.close(fig)
fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
axes[0].hist(nyears, bins=np.arange(0.5, len(years) + 1.5), color="tab:blue"); axes[0].set_xlabel("seasons in which the site held water"); axes[0].set_ylabel("sites")
axes[1].scatter(reg.area_km2, nyears + np.random.uniform(-0.2, 0.2, len(reg)), s=6, alpha=0.5); axes[1].set_xscale("log"); axes[1].set_xlabel("site area km2"); axes[1].set_ylabel("seasons with water")
bins = [0.05, 0.1, 0.2, 0.5, 1, 2, 5, 20]; cat = pd.cut(reg.area_km2, bins); fr = pd.Series(cent_in_sink).groupby(cat, observed=True).mean()
axes[2].bar(range(len(fr)), fr.values * 100, color="tab:green"); axes[2].set_xticks(range(len(fr))); axes[2].set_xticklabels([str(c) for c in fr.index], rotation=30, fontsize=8); axes[2].set_ylabel("% of sites with centroid in a DEM depression"); axes[2].set_ylim(0, 100)
fig.suptitle(f"Tile {TILE} sites: persistence, size, and the DEM correlation (reported, not enforced)"); fig.tight_layout(); fig.savefig(os.path.join(OUT, f"09_sites_{TILE}_stats.png"), dpi=100)
print("wrote", f"out/09_sites_{TILE}.{{geojson,txt}}, 09_registry, 09_site_years, 09_dunmire_crosswalk, and three figures")
