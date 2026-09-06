"""Figures from 12a: per window, the two chosen scenes side by side with the site outlines in red and the serial number of
each site (the G-number, without the coordinates). Writes out/12_truecolor/{TILE}_w{k}.png.
Run:  TILE=19_39 nice -n 15 $(cat .python_env) scripts/12b_truecolor_figs.py"""
import os, json, numpy as np, pandas as pd, geopandas as gpd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from shapely.geometry import box
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(ROOT, "out"); TC = os.path.join(OUT, "12_truecolor"); RAW = os.path.join(TC, "raw")
TILE = os.environ.get("TILE", "19_39")
g = gpd.read_file(os.path.join(OUT, f"09_sites_{TILE}.geojson")).to_crs(3413)
wins = json.load(open(os.path.join(RAW, f"{TILE}_windows.json")))
for w in wins:
    wx0, wy0, wx1, wy1 = w["window"]; sub = g[g.geometry.intersects(box(wx0, wy0, wx1, wy1))]
    seas = [y for y in w["seasons"] if os.path.exists(os.path.join(RAW, f"{TILE}_w{w['k']}_{y}.npy"))]
    fig, axes = plt.subplots(1, len(seas), figsize=(9 * len(seas), 9.3), squeeze=False); axes = axes[0]
    for ax, y in zip(axes, seas):
        ax.imshow(np.load(os.path.join(RAW, f"{TILE}_w{w['k']}_{y}.npy")), extent=(wx0, wx1, wy0, wy1))
        sub[~sub.ice_marginal].plot(ax=ax, facecolor="none", edgecolor="red", linewidth=1.0)
        if sub.ice_marginal.any(): sub[sub.ice_marginal].plot(ax=ax, facecolor="none", edgecolor="orange", linewidth=1.0, linestyle="--")
        for r in sub.itertuples():
          parts = list(r.geometry.geoms) if hasattr(r.geometry, "geoms") else [r.geometry]; big = max(q.area for q in parts)
          for part in [q for q in parts if q.area >= 0.2 * big]:  # label the substantial parts only
            c = part.representative_point(); ax.annotate(f"{r.serial}", (c.x, c.y), fontsize=7, color="yellow", ha="center", va="center", fontweight="bold",
                                                              bbox=dict(boxstyle="round,pad=0.15", fc="black", ec="none", alpha=0.55))
        ax.set_xlim(wx0, wx1); ax.set_ylim(wy0, wy1); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(f"{w['dates'][str(y)]}  (Sentinel-2 true colour)", fontsize=11)
    axes[0].plot([wx0 + 300, wx0 + 2300], [wy0 + 300, wy0 + 300], color="white", linewidth=3); axes[0].text(wx0 + 1300, wy0 + 450, "2 km", color="white", ha="center", fontsize=9)
    fig.suptitle(f"Tile {TILE}, window {w['k']} ({(wx1-wx0)/1000:.0f} km): red = site outline, number = serial of the Lake ID; {len(sub)} sites", fontsize=12, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(os.path.join(TC, f"{TILE}_w{w['k']}.png"), dpi=110); plt.close(fig)
    print("wrote", f"{TILE}_w{w['k']}.png", len(sub), "sites", seas)
