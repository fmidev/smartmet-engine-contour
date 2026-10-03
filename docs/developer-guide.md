# Contour engine developer guide

This guide is for developers who change `smartmet-engine-contour`, or call it from a
plugin. The engine turns a data matrix and its coordinates into isolines or isobands as
OGR geometries, with caching. The contouring algorithm itself is in the
[trax](https://github.com/fmidev/smartmet-library-trax) library (see its `docs/`).

[CLAUDE.md](../CLAUDE.md) has the architecture and cache details.

## Contents

1. [Building and testing](#1-building-and-testing)
2. [The API](#2-the-api)
3. [Options](#3-options)
4. [What contour() does](#4-what-contour-does)
5. [Caches](#5-caches)
6. [Configuration](#6-configuration)
7. [Compatibility](#7-compatibility)
8. [Known pitfalls](#8-known-pitfalls)

---

## 1. Building and testing

```bash
make
make test                               # tframe regression tests in test/
cd test && make GeosToolsTest && ./GeosToolsTest
```

The tests link the local `contour.so` and the installed querydata engine; `EngineTest`
loads both through `test/cnf/reactor.conf` with the querydata files in `test/data/`.

## 2. The API

```cpp
std::vector<OGRGeometryPtr> contour(std::size_t dataHash,
                                    const Fmi::SpatialReference& outputCRS,
                                    const NFmiDataMatrix<float>& values,
                                    const Fmi::CoordinateMatrix& coordinates,
                                    [const Fmi::Box& clipBox,]
                                    const Options& options) const;

std::vector<OGRGeometryPtr> crossection(NFmiFastQueryInfo& q, [const Spine::Parameter& z,]
                                        const Options& options,
                                        double lon1, double lat1, double lon2, double lat2,
                                        std::size_t steps) const;
```

* `contour()` returns one geometry per isovalue (isolines) or per range (isobands), in
  `outputCRS`. The coordinates are the grid points in that CRS (the querydata engine's
  `getWorldCoordinates*` gives them). A grid stored in the opposite row order is flipped
  internally.
* **`dataHash`** identifies the data: the caller computes it (WMS starts from
  `Engine::Querydata::hash_value(q)` and adds the options of computed fields). It is part
  of the cache key, so it **must** change whenever the values change.
* With a clip box (a map tile in output coordinates), only the cells overlapping the box are
  contoured.
* `crossection()` contours a vertical cross-section along the line, over the data's levels
  or, in the second form, with the `z` parameter as the vertical coordinate (used by the
  cross_section plugin).
* `getCacheSizes()`, `clearCache()` for admin reporting.

## 3. Options

`Options` is constructed from a list of isovalues (isolines) or of `Range`s
(`lolimit` / `hilimit`, open ends allowed) for isobands, and carries:

| Field | Effect |
|-------|--------|
| `interpolation` | Linear (default), midpoint or logarithmic interpolation along cell edges. |
| `extrapolation` | Grow the valid data into missing cells by this many steps, for data that covers only land and would otherwise leave gaps at the shoreline. |
| `multiplier`, `offset` | Unit conversion applied to the data first. |
| `smoother` (box, median, morphology) or `filter_size` / `filter_degree` (Savitzky-Golay) | Smooth the data before contouring; `smoother` wins if both are set. |
| `minarea` | Drop polygons smaller than this (km²). |
| `closed_range` | Whether the last isoband includes its upper limit. |
| `validate`, `strict`, `desliver` | OGC-validate the result; fail instead of repairing; remove slivers. |
| `subdivide` | Subdivide cells for smoother output. |
| `threads` | Contour in parallel (band-parallel in trax); defaults to `contour.threads`. |
| `parameter`, `time`, `level`, `bbox` | Descriptive; part of the options' hash. |

## 4. What contour() does

1. Look the result up in the contour cache (data hash + CRS + clip box + options).
2. Apply multiplier and offset, then smoothing, to get the grid to contour.
3. Choose a grid wrapper: `NormalGrid`; `PaddedGrid` (surrounded by missing values, so that
   open-ended isobands close at the grid edge); or `ShiftedGrid` (global data wrapping
   around the date line).
4. Get the **valid-cells mask** for the clip box from its cache (it depends only on the
   coordinates, the CRS and the box, so all parameters and times of a tile share it).
5. Run trax over the valid cells, build the geometries, and optionally validate,
   desliver and filter by area.
6. Cache and return the geometries.

## 5. Caches

| Cache | Key | Size |
|-------|-----|------|
| Contours | data hash + CRS + clip box + options | `cache.max_contours` entries |
| Coordinate analysis | coordinates + CRS | 1000 entries (fixed) |
| Valid cells | coordinates + CRS + clip box | `cache.max_valid_cells_mbytes` (256), over 16 shards |

The valid-cells cache matters for tiled requests: building a mask scans the whole grid, so
without the cache every tile costs as much as the full map. A mask larger than 1/16 of the
cache size is never cached.

## 6. Configuration

`cache.max_contours`, `cache.max_valid_cells_mbytes`, and `contour.threads` (default for
`Options::threads`).

## 7. Compatibility

`contour()` and `crossection()` are **non-virtual**, and plugins pass `Options` by
reference. Adding a field to `Options` changes its layout, so every plugin that builds
`Options` (wms, cross_section, …) must be rebuilt with the new engine.

## 8. Known pitfalls

* **The caller owns the data hash.** A hash that does not change when the data changes
  (for example one that omits the time or level) returns stale contours from the cache.
* **Infinite values** in the data produce NaN coordinates in trax; replace ±inf with
  missing values before contouring.
* **Large masks are not cached** (§5): very large grids with a small cache pay the full
  scan for every tile.
* **`Options` layout is ABI** (§7).
