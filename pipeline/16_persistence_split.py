"""Persistence surface + split rule (2026-09-07, Josh's "two neighbouring lakes that touch in a wet
season should not necessarily become one site").

The union of ten seasonal outlines is a binary mask: once two lakes touch in one wet season, no
geometry rule can recover that they were apart in the other nine. The information that can is the
per-pixel PERSISTENCE we already compute and then discard --- for each pixel, in how many of the ten
seasons was it inside an outline (0..10).

  nsea[y,x] = # seasons in which the pixel is inside a kept seasonal outline (same recipe as 08e)

On that surface both of Josh's failure modes are one mechanism:
  ice lid          the lid pixels are water in most seasons, occluded only in some scenes of one
                   season -> high persistence, no saddle, the patches stay one site.
  two lakes kissing the neck is water in 1-2 of 10 seasons while both cores are 8-10
                   -> a deep saddle -> two sites.

The test is a max-tree (component tree) prominence, the same device as MIN_PROM=3 on the DEM in
05_hybrid_ids_CW.py, but on observed water instead of topography --- so it also works for the 9-12 %
of lakes that sit in no ArcticDEM depression, and the DEM stays reported, never enforced (spec 3).

  threshold nsea at t = 10, 9, ... 1 and track components.  When two components merge at level t,
  the taller one continues and the other gets prominence = its peak level - t.  A site is more than
  one lake if >= 2 cores survive with prominence >= PROM and core size >= MIN_CORE_PX.

Separation of concerns: the appendage rule (08e) owns GAPS between disconnected union pieces; this
rule owns NECKS inside one connected piece.  So cores are counted within each connected component of
the site, and the scoring below reports single-piece sites (where this rule is the only decider)
apart from multi-piece ones (where the appendage rule already ruled).

Outputs
  out/16_nsea_{TILE}.npy            the persistence surface (cached; delete to rebuild, ~5 min)
  out/16_cores_{TILE}.csv           per site: pieces, cores, peaks, prominences, core sizes
  out/16_split_scores_{TILE}.txt    scored against Josh's and Claude's labels, PROM/MIN_CORE sweep
  out/16_split_cases_{TILE}.png     panels for the sites that split (or NAMED=... )

Run:  TILE=29_45 nice -n 15 $(cat .python_env) scripts/16_persistence_split.py
      PROM=3 MIN_CORE_PX=500 NAMED="G00212,G00113" ...
"""
import os, json, glob, time, numpy as np, pandas as pd, geopandas as gpd, scipy.ndimage as ndi
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from rasterio.features import rasterize
from affine import Affine
from shapely.geometry import box

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(ROOT, "out")
TILE = os.environ.get("TILE", "29_45")
PROM = int(os.environ.get("PROM", "3"))                    # seasons a core must stand above its saddle
MIN_CORE_PX = int(os.environ.get("MIN_CORE_PX", "500"))    # 0.05 km2 at 10 m, the outline threshold
MIN_PX = 500; FILL = 0.5; S8 = np.ones((3, 3), bool)
NAMED = [s for s in os.environ.get("NAMED", "").split(",") if s]

def disk(r):
    y, x = np.ogrid[-r:r + 1, -r:r + 1]; return (x * x + y * y) <= r * r

pref = f"08_s2counts_{TILE}_"
years = sorted(int(os.path.basename(f)[len(pref):].split("_")[0]) for f in glob.glob(os.path.join(OUT, pref + "*_meta.json")))
meta = json.load(open(os.path.join(OUT, f"{pref}{years[0]}_meta.json")))
tr = Affine(*meta["transform"]); N = meta["shape"][0]
x0, y0, x1, y1 = meta["x0"], meta["y0"], meta["x1"], meta["y1"]
lines = [f"Persistence split rule, tile {TILE}, seasons {years} — {time.strftime('%Y-%m-%d %H:%M %Z')}",
         f"PROM = {PROM} seasons, MIN_CORE_PX = {MIN_CORE_PX} ({MIN_CORE_PX * 1e-4:.2f} km2 at 10 m)"]
def say(s): print(s); lines.append(s)

# ---------------------------------------------------------------- 1. the persistence surface
nsf = os.path.join(OUT, f"16_nsea_{TILE}.npy")
if os.path.exists(nsf):
    nsea = np.load(nsf); say(f"persistence surface: cached {os.path.basename(nsf)}")
else:
    im = np.load(os.path.join(OUT, f"{TILE}_icemask_150m.npy")); imm = json.load(open(os.path.join(OUT, f"{TILE}_icemask_150m_meta.json")))
    xc = x0 + 10 * (np.arange(N) + 0.5); yc = y1 - 10 * (np.arange(N) + 0.5)
    ice = im[np.clip(np.rint((imm["y_row0"] - yc) / 150).astype(int), 0, im.shape[0] - 1)[:, None],
             np.clip(np.rint((xc - imm["x_col0"]) / 150).astype(int), 0, im.shape[1] - 1)[None, :]]
    nsea = np.zeros((N, N), np.uint8)
    for y in years:  # per-year outline recipe identical to 08e / 11a
        arr = np.where(ice, np.load(os.path.join(OUT, f"{pref}{y}_n_w50_toa.npy")), 0)
        w = ndi.binary_closing(arr >= 1, structure=disk(5)); lab, _ = ndi.label(w, structure=S8)
        cnt = np.bincount(lab.ravel()); ids = np.flatnonzero(cnt >= MIN_PX); ids = ids[ids > 0]
        keep = ids[ndi.mean((arr >= 1).astype(float), lab, ids) >= FILL]
        keep = keep[ndi.maximum(ndi.binary_erosion(np.isin(lab, keep), structure=disk(2)).astype(np.uint8), lab, keep).astype(bool)]
        nsea += np.isin(lab, keep).astype(np.uint8)
        say(f"  {y}: {len(keep)} outlines")
    np.save(nsf, nsea); say(f"persistence surface written {os.path.basename(nsf)}")
say(f"persistence: {(nsea > 0).sum() * 1e-4:.1f} km2 ever water; "
    + ", ".join(f"{t}+:{(nsea >= t).sum() * 1e-4:.0f}" for t in (1, 3, 5, 8, 10)) + " km2")

# ---------------------------------------------------------------- 2. sites, rasterised from the registry
sites = gpd.read_file(os.path.join(OUT, f"09_sites_{TILE}.geojson")).to_crs("EPSG:3413")  # registry is 4326; the grid is 3413
key = "site_id"   # geojson also has a "site" column, but that is a 1-based row index
sites = sites.reset_index(drop=True); sites["_i"] = np.arange(1, len(sites) + 1)
slab = rasterize(((g, i) for g, i in zip(sites.geometry, sites["_i"])), out_shape=(N, N), transform=tr, dtype=np.int32)
say(f"sites: {len(sites)} rasterised from 09_sites_{TILE}.geojson")

# ---------------------------------------------------------------- 3. max-tree prominence, per site
# The component tree does not depend on PROM/MIN_CORE_PX — those only filter its result — so it is
# built once per site and the sweep below just re-filters.  Bounding boxes come from find_objects so
# no step ever scans the full 10 000^2 raster per site.
def tree_of(v):
    """v = persistence values on one connected piece (0 elsewhere).  Returns
    (tallest_peak, tallest_size, absorbed) where absorbed = [(peak, saddle, prominence, core_px), ...]
    for every component that merged into a taller one."""
    top = int(v.max())
    if top == 0: return (0, 0, [])
    prev_ids = np.array([], int); prev_lab = None; prev_peak = {}; prev_size = {}; absorbed = []
    for t in range(top, 0, -1):
        lab, n = ndi.label(v >= t, structure=S8)
        if n == 0: continue
        sizes = np.bincount(lab.ravel(), minlength=n + 1)
        groups = {}
        if len(prev_ids):
            host = np.atleast_1d(ndi.maximum(lab, prev_lab, index=prev_ids)).astype(int)
            for pc, hc in zip(prev_ids, host): groups.setdefault(int(hc), []).append(int(pc))
        peak, size = {}, {}
        for c in range(1, n + 1):
            g = groups.get(c, [])
            if not g:                                             # a new peak appears at this level
                peak[c] = t; size[c] = int(sizes[c])
            else:
                g.sort(key=lambda pc: (prev_peak[pc], prev_size[pc]), reverse=True)
                peak[c] = prev_peak[g[0]]; size[c] = int(sizes[c])
                for pc in g[1:]:                                  # absorbed by a taller peak at level t
                    absorbed.append((prev_peak[pc], t, prev_peak[pc] - t, prev_size[pc]))
        prev_lab, prev_ids, prev_peak, prev_size = lab, np.arange(1, n + 1), peak, size
    tallp = max(prev_peak.values()) if prev_peak else top
    talls = max(prev_size.values()) if prev_size else 0
    return (tallp, talls, absorbed)

objs = ndi.find_objects(slab)
TREE = {}
for _, s_ in sites.iterrows():
    i = int(s_["_i"]); sl = objs[i - 1]
    if sl is None: continue
    m = (slab[sl] == i); v = np.where(m, nsea[sl], 0)
    pl, npieces = ndi.label(m, structure=S8)
    pieces = [tree_of(np.where(pl == p, v, 0)) for p in range(1, npieces + 1)]
    TREE[s_[key]] = dict(npieces=npieces, pieces=pieces, slice=sl)
say(f"component trees built for {len(TREE)} sites")

def cores_per_piece(sid, prom, min_core):
    """number of surviving cores in each connected piece of the site"""
    t = TREE.get(sid)
    if not t: return [0]
    return [1 + sum(1 for a in ab if a[2] >= prom and a[3] >= min_core) for (_, _, ab) in t["pieces"]]

rows = []
for sid, t in TREE.items():
    per = cores_per_piece(sid, PROM, MIN_CORE_PX)
    surv = [a for (_, _, ab) in t["pieces"] for a in ab if a[2] >= PROM and a[3] >= MIN_CORE_PX]
    rows.append(dict(site_id=sid, npieces=t["npieces"], n_cores=int(sum(per)),
                     max_cores_in_a_piece=int(max(per)), splits_a_piece=int(max(per) >= 2),
                     tall_peaks="|".join(str(tp) for (tp, _, _) in t["pieces"]),
                     split_peaks="|".join(str(a[0]) for a in surv),
                     split_saddles="|".join(str(a[1]) for a in surv),
                     split_proms="|".join(str(a[2]) for a in surv),
                     split_core_px="|".join(str(a[3]) for a in surv)))
cores = pd.DataFrame(rows)
cores.to_csv(os.path.join(OUT, f"16_cores_{TILE}.csv"), index=False)
say(f"cores: {len(cores)} sites; {int(cores.splits_a_piece.sum())} have >= 2 persistent cores inside one "
    f"connected piece at PROM={PROM} ({cores.splits_a_piece.mean() * 100:.1f} %)")

# ---------------------------------------------------------------- 4. score against the labels
# Both tiles were rebuilt (appendage rule) AFTER the closing test was labelled, so serials shifted and
# case_id no longer matches a current site_id.  Re-link each case to the current site by geometry:
# the site with the largest overlap with the case window (11a wrote it in EPSG:3413).
lf = os.path.join(OUT, f"11_labels_both_{TILE}.csv"); cf = os.path.join(OUT, f"11_cases_{TILE}.csv")
if os.path.exists(lf) and os.path.exists(cf):
    lab = pd.read_csv(lf).merge(pd.read_csv(cf)[["case_id", "xmin", "ymin", "xmax", "ymax"]], on="case_id", how="left")
    lab = lab[lab.pool.isin(["joined", "lid"])].dropna(subset=["xmin"])
    cw = gpd.GeoDataFrame(lab[["case_id"]].copy(),
                          geometry=[box(a, b, c, d) for a, b, c, d in zip(lab.xmin, lab.ymin, lab.xmax, lab.ymax)],
                          crs="EPSG:3413")
    ov = gpd.overlay(cw, sites[[key, "geometry"]], how="intersection", keep_geom_type=False)
    ov["ia"] = ov.area
    best = ov.sort_values("ia").groupby("case_id").tail(1)[["case_id", key]]
    n_before = len(lab)
    lab = lab.merge(best, on="case_id", how="inner").merge(cores, on=key, how="inner")
    L = lab[lab.josh.isin([1, 2])].copy()
    say(f"\nscoring set: {len(L)} of {n_before} joined/lid cases re-linked to a current site with a 1/2 label "
        f"(josh: {int((L.josh == 2).sum())} say two lakes, {int((L.josh == 1).sum())} say one)")
    say("\n  the rule fires only inside a connected piece; multi-piece sites were already decided by the")
    say("  appendage rule (gaps), so single-piece cases are where this rule is the only decider.")
    say("\n  PROM MIN_CORE |  single-piece n  acc_josh acc_claude |  all joined+lid n  acc_josh")
    best_s = None
    for prom in (1, 2, 3, 4, 5):
        for mc in (500, 1000, 2500):
            pred = {sid: max(cores_per_piece(sid, prom, mc)) >= 2 for sid in TREE}
            L["pred2"] = L[key].map(pred).fillna(False)
            S = L[L.npieces == 1]
            aj = (S.pred2 == (S.josh == 2)).mean() if len(S) else np.nan
            ac = (S.pred2 == (S.claude == 2)).mean() if len(S) else np.nan
            aa = (L.pred2 == (L.josh == 2)).mean() if len(L) else np.nan
            say(f"  {prom:4d} {mc:8d} | {len(S):15d} {aj:9.3f} {ac:10.3f} | {len(L):17d} {aa:9.3f}")
            if len(S) and (best_s is None or aj > best_s[0]): best_s = (aj, prom, mc)
    if best_s: say(f"\n  best single-piece accuracy vs Josh: {best_s[0]:.3f} at PROM={best_s[1]}, MIN_CORE_PX={best_s[2]}")
    # the baseline any split rule has to beat: never split anything
    Sb = L[L.npieces == 1]
    base = (Sb.josh == 1).mean() if len(Sb) else np.nan
    say(f"  BASELINE, never split: {base:.3f} vs Josh, {(Sb.claude == 1).mean():.3f} vs Claude "
        f"({int((Sb.josh == 2).sum())} of {len(Sb)} single-piece cases are really two lakes)")

    pred = {sid: max(cores_per_piece(sid, PROM, MIN_CORE_PX)) >= 2 for sid in TREE}
    L["pred2"] = L[key].map(pred).fillna(False); S = L[L.npieces == 1]
    tp = int((S.pred2 & (S.josh == 2)).sum()); fp = int((S.pred2 & (S.josh == 1)).sum())
    fn = int((~S.pred2 & (S.josh == 2)).sum()); tn = int((~S.pred2 & (S.josh == 1)).sum())
    say(f"\nat PROM={PROM}, MIN_CORE_PX={MIN_CORE_PX}, single-piece cases (n={len(S)}):")
    say(f"  correctly split {tp}   correctly kept one {tn}   over-split {fp}   missed splits {fn}")
    for nm, sub in (("missed", S[~S.pred2 & (S.josh == 2)]), ("over-split", S[S.pred2 & (S.josh == 1)])):
        if len(sub): say(f"  {nm}: " + ", ".join(f"{r[1].split('_')[0]}({r[2]})" for r in sub[[key, "pool"]].itertuples()))
    say("\n  flat-surface cases (the whole site exists in <= 2 seasons, so persistence cannot see a neck):")
    flat = S[S.tall_peaks.astype(str).str.split("|").str[0].astype(int) <= 2]
    say(f"    {len(flat)} of {len(S)}; josh calls {int((flat.josh == 2).sum())} of them two lakes")
    lab.to_csv(os.path.join(OUT, f"16_scored_{TILE}.csv"), index=False)
else:
    say(f"\n(no {os.path.basename(lf)} / {os.path.basename(cf)} — scoring skipped)")

# ---------------------------------------------------------------- 5. panels
all_ids = list(TREE)
show = [i for n in NAMED for i in all_ids if i.startswith(n)] if NAMED else cores[cores.splits_a_piece == 1].site_id.tolist()[:12]
if show:
    n = len(show); ncol = min(4, n); nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 3.9 * nrow), squeeze=False)
    cmap = plt.get_cmap("viridis", len(years) + 1); norm = BoundaryNorm(np.arange(-0.5, len(years) + 1.5), cmap.N)
    for ax, sid in zip(axes.ravel(), show):
        t = TREE[sid]; i = int(sites.loc[sites[key] == sid, "_i"].iloc[0]); sl = t["slice"]
        pad = 30
        r0 = max(0, sl[0].start - pad); r1 = min(N, sl[0].stop + pad)
        c0 = max(0, sl[1].start - pad); c1 = min(N, sl[1].stop + pad)
        m = slab[r0:r1, c0:c1] == i; v = np.where(m, nsea[r0:r1, c0:c1], 0)
        im_ = ax.imshow(v, cmap=cmap, norm=norm, interpolation="nearest")
        ax.contour(m.astype(float), levels=[0.5], colors="w", linewidths=0.6)
        per = cores_per_piece(sid, PROM, MIN_CORE_PX)
        surv = [a for (_, _, ab) in t["pieces"] for a in ab if a[2] >= PROM and a[3] >= MIN_CORE_PX]
        pk = "|".join(str(tp) for (tp, _, _) in t["pieces"])
        sp = ", ".join(f"peak {a[0]} over saddle {a[1]} (prom {a[2]}, {a[3] * 1e-4:.2f} km2)" for a in surv) or "none"
        ax.set_title(f"{sid.split('_')[0]}  pieces {t['npieces']}, cores {sum(per)}\ntall peak {pk}\nsecond core: {sp}", fontsize=6.5)
        ax.set_xticks([]); ax.set_yticks([])
    for ax in axes.ravel()[len(show):]: ax.axis("off")
    fig.colorbar(im_, ax=axes, shrink=0.6, label="seasons in which the pixel was inside an outline")
    f = os.path.join(OUT, f"16_split_cases_{TILE}.png"); fig.savefig(f, dpi=130, bbox_inches="tight"); plt.close(fig)
    say(f"\nfigure: {os.path.basename(f)} ({len(show)} panels)")

open(os.path.join(OUT, f"16_split_scores_{TILE}.txt"), "w").write("\n".join(lines) + "\n")
print(f"\nwrote out/16_split_scores_{TILE}.txt, out/16_cores_{TILE}.csv")
