%define DIRNAME contour
%define LIBNAME smartmet-%{DIRNAME}
%define SPECNAME smartmet-engine-%{DIRNAME}
Summary: SmartMet contour engine
Name: %{SPECNAME}
Version: 22.5.5
Release: 1%{?dist}.fmi
License: MIT
Group: SmartMet/Engines
URL: https://github.com/fmidev/smartmet-engine-contour
Source0: %{name}.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)
BuildRequires: boost169-devel
BuildRequires: bzip2-devel
BuildRequires: gcc-c++
%if %{defined el7}
BuildRequires: devtoolset-7-gcc-c++
#TestRequires: devtoolset-7-gcc-c++
%endif
BuildRequires: gdal34-devel
BuildRequires: geos310-devel
BuildRequires: make
BuildRequires: rpm-build
BuildRequires: libconfig17-devel
BuildRequires: smartmet-library-gis-devel >= 22.5.4
BuildRequires: smartmet-library-macgyver-devel >= 22.3.28
BuildRequires: smartmet-library-trax-devel >= 22.5.5
BuildRequires: smartmet-library-spine-devel >= 22.4.29
BuildRequires: sparsehash-devel
BuildRequires: zlib-devel
Requires: boost169-date-time
Requires: boost169-filesystem
Requires: boost169-iostreams
Requires: boost169-system
Requires: boost169-thread
Requires: gdal34-libs
Requires: geos310
Requires: smartmet-library-gis >= 22.5.4
Requires: smartmet-library-trax >= 22.5.5
Requires: smartmet-library-macgyver >= 22.3.28
Requires: smartmet-library-newbase >= 22.5.4
Requires: smartmet-library-spine >= 22.4.29
Requires: smartmet-library-timeseries >= 22.3.18
Requires: libconfig17

Provides: %{SPECNAME}
Obsoletes: smartmet-brainstorm-contour < 16.11.1
Obsoletes: smartmet-brainstorm-contour-debuginfo < 16.11.1
#TestRequires: boost169-devel
#TestRequires: bzip2-devel
#TestRequires: gcc-c++
#TestRequires: gdal34-devel
#TestRequires: geos310-devel
#TestRequires: libjpeg-turbo-devel
#TestRequires: libpng-devel
#TestRequires: smartmet-engine-querydata >= 21.1.21
#TestRequires: smartmet-engine-querydata-devel >= 21.1.21
#TestRequires: smartmet-library-regression >= 21.1.21
#TestRequires: smartmet-library-spine-devel >= 22.4.29
#TestRequires: smartmet-library-timeseries-devel >= 22.3.18
#TestRequires: smartmet-library-trax >= 22.5.5
#TestRequires: smartmet-library-trax-devel >= 22.5.5
#TestRequires: smartmet-test-data >= 20.10.29
#TestRequires: zlib-devel


%description
SmartMet contour engine

%package -n %{SPECNAME}-devel
Summary: SmartMet %{SPECNAME} development headers
Group: SmartMet/Development
Provides: %{SPECNAME}-devel
Requires: %{SPECNAME} = %{version}-%{release}
Obsoletes: smartmet-brainstorm-contour-devel < 16.11.1
%description -n %{SPECNAME}-devel
SmartMet %{SPECNAME} development headers.

%prep
rm -rf $RPM_BUILD_ROOT

%setup -q -n %{SPECNAME}
 
%build -q -n %{SPECNAME}
make %{_smp_mflags}

%install
%makeinstall
mkdir -p $RPM_BUILD_ROOT%{_sysconfdir}/smartmet/engines

%clean
rm -rf $RPM_BUILD_ROOT

%files -n %{SPECNAME}
%defattr(0775,root,root,0775)
%{_datadir}/smartmet/engines/%{DIRNAME}.so

%files -n %{SPECNAME}-devel
%defattr(0664,root,root,0775)
%{_includedir}/smartmet/engines/%{DIRNAME}

%changelog
* Thu May  5 2022 Mika Heiskanen <mika.heiskanen@fmi.fi> - 22.5.5-1.fmi
- Repackaged since Trax ABI switched from using doubles to floats

* Wed May  4 2022 Mika Heiskanen <mika.heiskanen@fmi.fi> - 22.5.4-2.fmi
- Use Trax for contouring
- Use tiled contouring by default
- Fixed data smoothing

* Fri Mar 11 2022 Mika Heiskanen <mika.heiskanen@fmi.fi> - 22.3.11-1.fmi
- Repackaged due to refactored library dependencies

* Fri Jan 21 2022 Andris Pavēnis <andris.pavenis@fmi.fi> 22.1.21-1.fmi
- Repackage due to upgrade of packages from PGDG repo: gdal-3.4, geos-3.10, proj-8.2

* Tue Dec  7 2021 Andris Pavēnis <andris.pavenis@fmi.fi> 21.12.7-1.fmi
- Update to postgresql 13 and gdal 3.3

* Tue Sep 28 2021 Andris Pavēnis <andris.pavenis@fmi.fi> 21.9.28-1.fmi
- Repackage due to dependency change: moving libconfig files to differentr directory

* Mon Sep 13 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.9.13-1.fmi
- Repackaged due to Fmi::Cache statistics fixes

* Tue Sep  7 2021 Andris Pavēnis <andris.pavenis@fmi.fi> 21.9.7-1.fmi
- Rebuild due to dependency changes

* Mon Aug 30 2021 Anssi Reponen <anssi.reponen@fmi.fi> - 21.8.30-1.fmi
- Cache counters added (BRAINSTORM-1005)

* Tue Aug 17 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.8.17-1.fmi
- Repackaged due to shutdown ABI changes

* Wed Jul 28 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.7.28-1.fmi
- Silenced compiler warnings

* Mon Jun 28 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.6.28-1.fmi
- Detect newbase version from external WGS84 define

* Thu May 20 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.5.20-2.fmi
- Repackaged with improved hashing functions

* Thu May 20 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.5.20-1.fmi
- Use Fmi hash functions, boost::hash_combine produces too many collisions

* Thu May  6 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.5.6-1.fmi
- Repackaged due to ABI changes in NFmiAzimuthalArea

* Tue Apr 13 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.4.13-1.fmi
- Repackaged with the latest Tron library: match isolines with isoband boundaries

* Fri Apr  9 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.4.9-3.fmi
- Avoid unnecessary copy of coordinates if no flipping is needed

* Fri Apr  9 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.4.9-2.fmi
- More Tron optimizations

* Thu Apr  8 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.4.8-1.fmi
- Rebuilt with optimized Tron library

* Tue Mar 23 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.3.23-1.fmi
- Repackaged due to geos39 updates

* Thu Feb 18 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.2.18-1.fmi
- Use new NFmiPoint API

* Thu Feb 11 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.2.11-1.fmi
- Use GIS-library definition of OGRGeometryPtr for consistency

* Wed Feb 10 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.2.10-1.fmi
- Use CoordinateMatrix APIs

* Thu Jan 14 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.1.14-1.fmi
- Repackaged smartmet to resolve debuginfo issues

* Tue Jan  5 2021 Mika Heiskanen <mika.heiskanen@fmi.fi> - 21.1.5-1.fmi
- Upgrade to GEOS 3.9

* Mon Dec 28 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.12.28-1.fmi
- Upgraded tron library

* Tue Dec 15 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.12.15-1.fmi
- Upgrade to pgdg12

* Fri Nov  6 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.11.6-1.fmi
- Modified includes to enable gdal 3.0 support

* Tue Oct  6 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.10.6-1.fmi
- Enable sensible relative libconfig include paths

* Wed Sep 23 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.8.23-1.fmi
- Use Fmi::Exception instead of Spine::Exception

* Thu Aug 27 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.8.27-1.fmi
- NFmiGrid API changed

* Wed Aug 26 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.8.26-1.fmi
- Numerous newbase API changes

* Fri Aug 21 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.8.21-1.fmi
- Upgraded to fmt 6.2

* Thu Jul  2 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.7.2-1.fmi
- SpatialReference API changed

* Wed May 13 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.5.13-1.fmi
- Repackaged since Spine Parameter class ABI changed

* Thu Apr 23 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.4.22-1.fmi
- Rewrote global grid wrapping code

* Wed Apr 22 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.4.20-1.fmi
- Improved gdal30/geos38 detection

* Sat Apr 18 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.4.18-1.fmi
- Upgraded to Boost 1.69

* Thu Apr  9 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.4.9-1.fmi
- Repackaged since code inlines from Tron changed

* Thu Apr  2 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.4.2-1.fmi
- Use NFmiCoordinateMatrix instead of NFmiDataMatrix<NFmiPoint>

* Mon Mar 30 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.3.30-1.fmi
- Full repackaging of the server

* Thu Feb 13 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.2.13-1.fmi
- Use GDAL forward declarations in Engine API to avoid dependency escalation

* Wed Feb 12 2020 Mika Heiskanen <mika.heiskanen@fmi.fi> - 20.2.12-1.fmi
- Repackaged due to API changes

* Fri Dec 13 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.12.13-1.fmi
- Repackaged due to NFmiArea API changes

* Wed Dec 11 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.12.11-1.fmi
- Upgrade to GDAL 3.0 and GEOS 3.8

* Wed Nov 20 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.11.20-1.fmi
- Repackaged due to newbase API changes

* Thu Oct 31 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.10.31-2.fmi
- Rebuild due to Tron API changes (shared_ptr)

* Thu Oct 31 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.10.31-1.fmi
- Rebuilt due to newbase API/ABI changes

* Thu Sep 26 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.9.26-1.fmi
- Added support for GEOS 3

* Fri Feb 22 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.2.22-1.fmi
- Fixed contour cache hash value to include the minarea setting

* Thu Feb 21 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.2.21-1.fmi
- Added optional minarea setting for closed contours

* Tue Feb 12 2019 Mika Heiskanen <mika.heiskanen@fmi.fi> - 19.2.12-1.fmi
- Added extrapolation support

* Wed Jul 25 2018 Mika Heiskanen <mika.heiskanen@fmi.fi> - 18.7.25-1.fmi
- Prefer nullptr over NULL

* Fri Jun 29 2018 Mika Heiskanen <mika.heiskanen@fmi.fi> - 18.6.29-1.fmi
- Use less threads for contouring to avoid high CPU load peaks

* Wed May  9 2018 Mika Heiskanen <mika.heiskanen@fmi.fi> - 18.4.9-1.fmi
- Optimized missing value handling for speed

* Sat Apr  7 2018 Mika Heiskanen <mika.heiskanen@fmi.fi> - 18.4.7-1.fmi
- Upgrade to boost 1.66

* Tue Mar 20 2018 Mika Heiskanen <mika.heiskanen@fmi.fi> - 18.3.20-1.fmi
- Full repackaging of the server

* Fri Feb  9 2018 Mika Heiskanen <mika.heiskanen@fmi.fi> - 18.2.9-1.fmi
- Repackaged since base class SmartMetEngine size changed

* Wed Nov  1 2017 Mika Heiskanen <mika.heiskanen@fmi.fi> - 17.11.1-1.fmi
- Rebuilt due to GIS-library API change

* Wed Sep 20 2017 Mika Heiskanen <mika.heiskanen@fmi.fi> - 17.9.20-1.fmi
- Conturing now uses -Inf and Inf to represent open intervals
- Disabled -Ofast

* Tue Sep 12 2017 Mika Heiskanen <mika.heiskanen@fmi.fi> - 17.9.12-1.fmi
- Added automatic detection of inverted grids to enforce geometry winding rules

* Sun Aug 27 2017 Mika Heiskanen <mika.heiskanen@fmi.fi> - 17.8.28-1.fmi
- Upgrade to boost 1.65
- Detect when missing level values when calculating cross-sections

* Tue May 30 2017 Mika Heiskanen <mika.heiskanen@fmi.fi> - 17.5.30-1.fmi
- Do not use maximal concurrency for conturing until better load balancing is implemented

* Tue Apr 25 2017 Mika Heiskanen <mika.heiskanen@fmi.fi> - 17.4.25-2.fmi
- Bug fix: return empty vector if the input data is zero size, not a vector of empty shared pointers

* Tue Apr 25 2017 Mika Heiskanen <mika.heiskanen@fmi.fi> - 17.4.25-1.fmi
- Prevent exceptions during contouring cause termination via thread exit

* Wed Apr 19 2017 Mika Heiskanen <mika.heiskanen@fmi.fi> - 17.4.19-1.fmi
- Contouring is now done in parallel
- Added Engine::clearCache to facilitate speed comparisons

* Mon Apr 10 2017 Mika Heiskanen <mika.heiskanen@fmi.fi> - 17.4.10-1.fmi
- Added boolean to the API which indicates whether the data needs to be wrapped around to cover the globe

* Wed Mar 15 2017 Mika Heiskanen <mika.heiskanen@fmi.fi> - 17.3.15-1.fmi
- Recompiled since Spine::Exception changed

* Sat Feb 11 2017 Mika Heiskanen <mika.heiskanen@fmi.fi> - 17.2.11-1.fmi
- Recpackaged due to newbase API change in function CreateNewArea

* Tue Feb  7 2017 Mika Heiskanen <mika.heiskanen@fmi.fi> - 17.2.7-1.fmi
- Tron fixes to isoline building algorithm

* Tue Jan 31 2017 Mika Heiskanen <mika.heiskanen@fmi.fi> - 17.1.31-1.fmi
- Recompiled with latest Tron without Google dense_hash for much better speed

* Mon Jan 30 2017 Mika Heiskanen <mika.heiskanen@fmi.fi> - 17.1.30-1.fmi
- Tron now handles self-touching isolines correctly

* Thu Jan 19 2017 Mika Heiskanen <mika.heiskanen@fmi.fi> - 17.1.19-1.fmi
- Recompiled with more robost isoline calculations

* Wed Jan 18 2017 Mika Heiskanen <mika.heiskanen@fmi.fi> - 17.1.18-1.fmi
- Bug fix to Tron in handling self touching isolines

* Wed Jan  4 2017 Mika Heiskanen <mika.heiskanen@fmi.fi> - 17.1.4-1.fmi
- Changed to use renamed SmartMet base libraries

* Wed Nov 30 2016 Mika Heiskanen <mika.heiskanen@fmi.fi> - 16.11.30-1.fmi
- New releace with refactored configuration files

* Fri Nov 18 2016 Mika Heiskanen <mika.heiskanen@fmi.fi> - 16.11.18-1.fmi
- Rebuild with extra safety checks during contouring

* Tue Nov  1 2016 Mika Heiskanen <mika.heiskanen@fmi.fi> - 16.11.1-1.fmi
- Namespace changed. Pimple class renamed to Impl.

* Tue Sep 20 2016 Mika Heiskanen <mika.heiskanen@fmi.fi> - 16.9.20-1.fmi
- Geometries inserted into cache one by one, not as a vector

* Tue Sep 13 2016 Mika Heiskanen <mika.heiskanen@fmi.fi> - 16.9.13-1.fmi
- Contour-engine API change: vector of isolines/isobands queried at once
- Hints-cache removed

* Tue Sep  6 2016 Mika Heiskanen <mika.heiskanen@fmi.fi> - 16.9.6-1.fmi
- New exception handler

* Mon Aug 15 2016 Mika Heiskanen <mika.heiskanen@fmi.fi> - 16.8.15-1.fmi
- Full recompile

* Tue Jun 14 2016 Mika Heiskanen <mika.heiskanen@fmi.fi> - 16.6.14-1.fmi
- Full recompile

* Thu Jun  2 2016 Mika Heiskanen <mika.heiskanen@fmi.fi> - 16.6.2-1.fmi
- Full recompile

* Wed Jun  1 2016 Mika Heiskanen <mika.heiskanen@fmi.fi> - 16.6.1-1.fmi
- Added graceful shutdown

* Wed May 11 2016 Mika Heiskanen <mika.heiskanen@fmi.fi> - 16.5.11-1.fmi
- Reduced cache size for data hints from 5000 to 1000 to reduce memory use
- Added cache size reporting

* Mon Apr 25 2016 Mika Heiskanen <mika.heiskanen@fmi.fi> - 16.4.25-1.fmi
- Improved isoline contouring

* Wed Apr 20 2016 Tuomo Lauri <tuomo.lauri@fmi.fi> - 16.4.20-1.fmi
- Removed QEngine dependencies

* Fri Apr  1 2016 Mika Heiskanen <mika.heiskanen@fmi.fi> - 16.4.1-1.fmi
- Using -Ofast for maximum speed
- Recompiled with Tron-library with a faster hole assignment algorithm

* Tue Mar 29 2016 Mika Heiskanen <mika.heiskanen@fmi.fi> - 16.3.29-1.fmi
- Recompiled with a faster Tron-contouring library

* Mon Mar  7 2016 Tuomo Lauri <tuomo.lauri@fmi.fi> - 16.3.4-1.fmi
- Compiled against new QEngine

* Mon Jan 18 2016 Mika Heiskanen <mika.heiskanen@fmi.fi> - 16.1.18-1.fmi
- newbase API changed, full recompile

* Mon Oct 26 2015 Mika Heiskanen <mika.heiskanen@fmi.fi> - 15.10.26-1.fmi
- Added proper debuginfo packaging

* Thu Oct 15 2015 Mika Heiskanen <mika.heiskanen@fmi.fi> - 15.10.15-1.fmi
- Recompile with the latest Tron-library with polygon builder fixes.

* Fri Sep 18 2015 Mika Heiskanen <mika.heiskanen@fmi.fi> - 15.9.18-1.fmi
- Handle pole and 180th meridian problems when reprojecting data

* Thu Sep 17 2015 Mika Heiskanen <mika.heiskanen@fmi.fi> - 15.9.17-1.fmi
- Flip cross-section grids if necessary to force a right handed coordinate system

* Wed Sep 16 2015 Mika Heiskanen <mika.heiskanen@fmi.fi> - 15.9.16-2.fmi
- Fixes to Tron rounding errors

* Wed Sep 16 2015 Mika Heiskanen <mika.heiskanen@fmi.fi> - 15.9.16-1.fmi
- Fixes to Tron saddle point handling

* Tue Sep 15 2015 Tuomo Lauri <tuomo.lauri@fmi.fi> - 15.9.15-1.fmi
- Contouring method now throws error if parameter is unavailable

* Mon Sep 14 2015 Mika Heiskanen <mika.heiskanen@fmi.fi> - 15.9.14-1.fmi
- Switch to using the faster FmiBuilder, GEOS builds polygons too slowly

* Tue Aug 18 2015 Mika Heiskanen <mika.heiskanen@fmi.fi> - 15.8.18-1.fmi
- Recompile forced by HTTP API changes

* Mon Aug 17 2015 Mika Heiskanen <mika.heiskanen@fmi.fi> - 15.8.17-1.fmi
- Use -fno-omit-frame-pointer to improve perf use

* Thu Apr  9 2015 Mika Heiskanen <mika.heiskanen@fmi.fi> - 15.4.9-1.fmi
- newbase API changed

* Wed Apr  8 2015 Mika Heiskanen <mika.heiskanen@fmi.fi> - 15.4.8-1.fmi
- Dynamic linking of smartmet libraries into use
- Added caching of querydata projections
- Added caching of querydata values
- Added caching of querydata value quadtrees

* Tue Feb 24 2015 Mika Heiskanen <mika.heiskanen@fmi.fi> - 15.2.24-1.fmi
- Recompiled with the latest newbase with symbol changes

* Wed Dec 17 2014 Mika Heiskanen <mika.heiskanen@fmi.fi> - 14.12.17-1.fmi
- Recompiled due to spine API changes

* Thu Nov 13 2014 Mika Heiskanen <mika.heiskanen@fmi.fi> - 14.11.13-1.fmi
- Recompiled due to newbase API changes

* Fri Oct 24 2014 Tuomo Lauri <tuomo.lauri@fmi.fi> - 14.10.24-1.fmi
- Decreased cache size to limit memory consumption

* Tue Sep  2 2014 Mika Heiskanen <mika.heiskanen@fmi.fi> - 14.9.2-1.fmi
- Fixed zparameter to work when it is the first parameter in the data (index 0)

* Fri Aug 29 2014 Mika Heiskanen <mika.heiskanen@fmi.fi> - 14.8.29-1.fmi
- Added caching of calclated contours

* Tue Aug 26 2014 Mika Heiskanen <mika.heiskanen@fmi.fi> - 14.8.26-1.fmi
- Added possibility take cross section with any parameter as a height value

* Thu Aug 21 2014 Tuomo Lauri <tuomo.lauri@fmi.fi> - 14.8.21-2.fmi
- Fixed spatialref reference-counting related memory leak

* Thu Aug 21 2014 Tuomo Lauri <tuomo.lauri@fmi.fi> - 14.8.21-1.fmi
- Compiled with the fixed Tron

* Thu Aug  7 2014 Mika Heiskanen <mika.heiskanen@fmi.fi> - 14.8.7-1.fmi
- Fixed a dereference of a null shared pointer

* Wed Aug  6 2014 Mika Heiskanen <mika.heiskanen@fmi.fi> - 14.8.6-2.fmi
- Fixed the contourer to use geographic coordinates for latlon and rotated latlon data

* Wed Jul 30 2014 Mika Heiskanen <mika.heiskanen@fmi.fi> - 14.7.30-1.fmi
- Added support for contouring cross sections

* Fri Jul 25 2014 Mika Heiskanen <mika.heiskanen@fmi.fi> - 14.7.25-1.fmi
- Rebuilt with the latest GDAL

* Wed Jul  2 2014 Mika Heiskanen <mika.heiskanen@fmi.fi> - 14.7.2-1.fmi
- Isobands can now be requested with boost::optional lolimit and hilimit

* Thu Jun 26 2014 Mika Heiskanen <mika.heiskanen@fmi.fi> - 14.6.16-1.fmi
- Added support for unit conversions
- Added support for selecting the level

* Wed Jun 11 2014 Mika Heiskanen <mika.heiskanen@fmi.fi> - 14.6.11-1.fmi
- Recompiled with latest libgis API

* Wed Jun  4 2014 Mika Heiskanen <mika.heiskanen@fmi.fi> - 14.6.4-1.fmi
- Refactored to use libgis

* Thu Feb 20 2014 Mika Heiskanen <mika.heiskanen@fmi.fi> - 14.2.20-1.fmi
- First alpha version doing actual contouring

* Wed Feb 12 2014 Mika Heiskanen <mika.heiskanen@fmi.fi> - 14.2.12-1.fmi
- First alpha version
