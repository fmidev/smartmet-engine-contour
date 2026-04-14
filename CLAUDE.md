# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

The contour engine (`smartmet-engine-contour`) generates isolines and isobands from gridded weather data on demand. It is loaded as a shared library by SmartMet Server and used primarily by the WMS plugin to render contoured weather maps.

## Build commands

```bash
make                # Build contour.so
make test           # Run all tests (cd test && make test)
make format         # clang-format all source and test files
make clean          # Remove build artifacts
make install        # Install headers to $(includedir)/smartmet/engines/contour/ and .so to $(enginedir)
make rpm            # Build RPM package
```

To run a single test:

```bash
cd test && make GeosToolsTest && ./GeosToolsTest
```

Tests compile with `-O0 -g` and link against the locally built `contour.so` plus the system-installed `querydata.so` engine. The `EngineTest` requires both engines loaded via a reactor configuration (`test/cnf/reactor.conf`).

## Architecture

**Namespace**: `SmartMet::Engine::Contour`

**Engine** (`Engine.h/.cpp`) — Public interface, inherits `Spine::SmartMetEngine`. Uses pimpl (`Engine::Impl`) to hide caches and internals. Dynamically loaded by the server via the C extern `engine_class_creator()` function. Key methods:
- `contour()` — produces `std::vector<OGRGeometryPtr>` isolines/isobands from a data matrix + coordinate matrix, with optional clip box
- `crossection()` — produces cross-section contours along a geographic line

**Options** (`Options.h/.cpp`) — All contouring parameters. Constructed with either `std::vector<double>` (isovalues for isolines) or `std::vector<Range>` (lo/hi limits for isobands). Includes data transformation (multiplier/offset), Savitzky-Golay smoothing (filter_size/filter_degree), minimum area filtering, OGC validation, and desliver options.

**Grid hierarchy** — Four grid implementations, all inheriting `BaseGrid` (which inherits `Trax::Grid`):
- `NormalGrid` — regular matrix with optional world-data wraparound
- `PaddedGrid` — surrounded by NaN-filled cells for missing-data isobands outside grid bounds
- `ShiftedGrid` — handles globe data with column wraparound shift
- `MirrorGrid` — extends grid at borders by mirroring trends (used by Savitzky-Golay filter)

**SavitzkyGolay2D** — 2D smoothing filter applied to grids before contouring. Size range 0-6, degree range 0-5.

**GeosTools** — Utility for converting GEOS geometries to SVG path format with precision control.

**Caching** (`Engine::Impl`) — Dual LRU cache: contour geometry cache and coordinate analysis cache. Hash-based keys. Cache size configured via libconfig (default 10,000 entries).

## Key dependencies

- **trax** — marching-squares contouring algorithm (the actual isoline/isoband computation)
- **gis** — coordinate matrices, spatial references, bounding boxes
- **spine** — engine base class, reactor, HTTP/config framework
- **newbase** — `NFmiDataMatrix<float>` (the gridded data format), `NFmiFastQueryInfo`
- **GEOS** — geometry operations (built with `-DUSE_UNSTABLE_GEOS_CPP_API`)
- **GDAL/Proj** — geospatial transformations

## Testing

Tests use the SmartMet `tframe.h` regression framework (not Boost.Test). Assertions use `TEST_FAILED()` macro. Test data lives in `test/data/` (SQD files) with config in `test/cnf/`.
