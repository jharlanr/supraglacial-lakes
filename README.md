# supraglacial-lakes

A living database of supraglacial lakes with persistent IDs. Greenland first; the design is not Greenland-specific.

**Status (2026-09-05):** design v0.0 agreed, recipe tested on one 100 km tile (ArcticDEM tile 19_39, central-west Greenland)
for two seasons; the ten-season build of that tile is in progress. Private until the first data release.

## What a lake is here

- **Outline (per year):** the seasonal maximum extent of water in one melt season — every Sentinel-2 scene June–September,
  10 m, water where the blue–red NDWI exceeds 0.5 (Dunmire et al. 2021, from Miles et al. 2017), no scene-level cloud filter,
  a small closing, connected components ≥ 0.05 km², and a fill filter that drops slush and crevasse fields.
- **Lake (site):** the union of the per-year outlines over a basis of seasons (v0.0: 2016–2025), closed with a ~150 m kernel.
  The union is the lake's footprint polygon; the DEM is reported per lake as an attribute, never enforced.
- **Observation:** water inside the footprint per scene date, plus spill beyond it and a shared flag when one body spans two lakes.

## The ID

`G00000_N691234_W0485678` — `G` for Greenland, a five-digit serial for quick referral, then the registration centroid
(69.1234 N, 048.5678 W) as hemisphere letter + fixed-width digits with four implied decimals. No signs or points.

Rules: a serial is never reused or reassigned; the first registry numbers lakes north to south; each later batch of
additions is appended after the last serial, north to south within the batch; the coordinate part is a label frozen at
registration and is never recomputed when a polygon is redrawn. New lakes in later (or earlier) years are registered by
the append rule; lakes are never merged or renumbered by it — merges happen only in a versioned rebuild with a crosswalk.

## Layout

- `docs/spec_v0.0.md` — the design, with the evidence behind each choice.
- `docs/literature_lake_definition_and_identity.md` — how 30 papers define a lake and "the same lake" across years.
- `pipeline/` — the scripts: `08a` submits an Earth Engine export of per-pixel water counts for one tile-season;
  `08b` downloads it; `08c` sweeps the outline knobs against Dunmire's outlines; `08d` tests multi-season sites;
  `08e` builds sites, the registry, DEM attributes, the Dunmire crosswalk and figures. They currently expect the
  exploration workspace's layout (`out/`, `../labels/dunmire/`, a `.gee_project` file); a `data/` layout for this repo is the next step.
- `env/gee.yml` — the conda env (earthengine-api, geopandas, rasterio).

## Precedents

Sites-from-union follows Selmes et al. 2011; sinks-with-IDs follows Hochreuther et al. 2021 and Lutz et al. 2025;
the water recipe follows Dunmire et al. 2021; the ID and versioning pattern follows GLIMS / the Randolph Glacier Inventory.
Licence for code and data: to be decided before the first public release.
