# smartmet-engine-contour

Part of [SmartMet Server](https://github.com/fmidev/smartmet-server). See the [SmartMet Server documentation](https://github.com/fmidev/smartmet-server) for a full overview of the ecosystem.

## Overview

The contour engine generates isolines and isobands from gridded weather data. It is used by the WMS plugin to render contoured weather maps on demand.

## Features

- Isoline (contour line) generation from 2D gridded scalar fields
- Isoband (filled contour) generation
- Caching of computed contours for performance
- Integration with [smartmet-library-trax](https://github.com/fmidev/smartmet-library-trax) for marching-squares contouring
- Optional bilinear cell subdivision (`Options.subdivide`, `int` 0..10) that replaces the straight marching-squares segment between cell edge intersections with `n-1` samples on the true bilinear level curve. Forwarded to `Trax::Contour::subdivide()` and included in the contour cache key, so different subdivide values produce distinct cached entries.

## License

MIT — see [LICENSE](LICENSE)

## Contributing

Bug reports and pull requests are welcome on [GitHub](../../issues).
