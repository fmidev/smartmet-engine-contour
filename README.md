# smartmet-engine-contour

Part of [SmartMet Server](https://github.com/fmidev/smartmet-server). See the [SmartMet Server documentation](https://github.com/fmidev/smartmet-server) for a full overview of the ecosystem.

## Overview

The contour engine generates isolines and isobands from gridded weather data. It is used by the WMS plugin to render contoured weather maps on demand.

## Features

- Isoline (contour line) generation from 2D gridded scalar fields
- Isoband (filled contour) generation
- Caching of computed contours for performance
- Integration with [smartmet-library-trax](https://github.com/fmidev/smartmet-library-trax) for marching-squares contouring

## License

MIT — see [LICENSE](LICENSE)

## Contributing

Bug reports and pull requests are welcome on [GitHub](../../issues).
