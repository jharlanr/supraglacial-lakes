# supraglacial-lakes

A living database of supraglacial lakes with persistent IDs. Greenland first; the design is not Greenland-specific.

**Status (2026-09-05):** design v0.0 agreed, recipe tested on one 100 km tile (ArcticDEM tile 19_39, central-west Greenland)
for two seasons; the ten-season build of that tile is in progress. Private until the first data release.

## What a lake is here

- **Outline (per year):** the seasonal maximum extent of water in one melt season — every Sentinel-2 scene June–September,
  10 m, water where the blue–red NDWI exceeds 0.5 (Dunmire et al. 2021, from Miles et al. 2017), no scene-level cloud filter,
  a small closing, connected components ≥ 0.05 km², a fill filter that drops slush and crevasse fields, a 50 m core
  (drops swath-edge lines), and a minimum width of 40 m (channels and pond-chain necks are cut, lake shapes untouched).
- **Lake (site):** the union of the per-year outlines over a basis of seasons (v0.0: 2016–2025), on the ice sheet
  (BedMachine mask; water that touches rock or ocean is an ice-marginal lake and is excluded), with separate bodies of water joined by the *appendage rule*: two bodies are one lake if they are
  within 50 m of each other (an ice lid splitting a lake), or within 250 m when the smaller is under one fifth of the
  larger (a tail or pond hanging off a lake, which does not merit its own ID). Otherwise they are separate lakes. The
  three numbers were set by a labelled test of 105 cases (`tests/closing_19_39/`). The union is the lake's water-seen
  footprint polygon (it under-draws lakes under a persistent ice lid); the DEM is reported per lake as an attribute,
  never enforced, along with review flags (`ice_marginal`, `thin`).
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
  `08e` builds sites (appendage rule), the registry, DEM attributes and review flags, the Dunmire crosswalk and figures; `08f` drives all ten seasons of one tile; `08g` cuts the ice-sheet domain mask from BedMachine v6; `11e` scores the appendage rule against the labels; `12a`/`12b` make the true-colour check figures. They currently expect the
  exploration workspace's layout (`out/`, `../labels/dunmire/`, a `.gee_project` file); a `data/` layout for this repo is the next step.
- `registry/v0.0-test/tile_19_39/` — the first test registry (one 100 km ArcticDEM tile, central-west Greenland,
  ten seasons 2016–2025): `09_registry_19_39.csv` (one row per Lake ID: serial, centroid, area, seasons present,
  DEM attributes), `09_sites_19_39.geojson` (the site polygons), `09_site_years_19_39.csv` (per-season presence
  and area), `09_dunmire_crosswalk_19_39.csv` (which Dunmire 2018/2019 lakes each site holds), `09_sites_19_39.txt`
  (the summary). Test output, not a release: IDs here will be reissued when the Greenland-wide v0.0 registry is built.
- `registry/v0.0-test/tile_29_45/` — the second test tile (NE Greenland, 78–79° N), the independence check.
- `tests/closing_19_39/` — the closing-accuracy test: 105 cases, both labellers' answers, the bridged-gap table, and
  `11_rule_scores_19_39.txt`, the scoring that fixed the appendage rule (50 m always; 250 m if the smaller body is under
  one fifth of the larger).
- `docs/figures/` — per tile: map, showcase and statistics figures, and `truecolor/` — six 12 km windows per tile with
  the site outlines and serials over Sentinel-2 true colour on two summers.
- `env/gee.yml` — the conda env (earthengine-api, geopandas, rasterio).

## Precedents

Sites-from-union follows Selmes et al. 2011; sinks-with-IDs follows Hochreuther et al. 2021 and Lutz et al. 2025;
the water recipe follows Dunmire et al. 2021; the ID and versioning pattern follows GLIMS / the Randolph Glacier Inventory.
Licence for code and data: to be decided before the first public release.
