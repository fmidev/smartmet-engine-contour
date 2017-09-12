%define DIRNAME contour
%define LIBNAME smartmet-%{DIRNAME}
%define SPECNAME smartmet-engine-%{DIRNAME}
Summary: SmartMet contour engine
Name: %{SPECNAME}
Version: 17.9.12
Release: 1%{?dist}.fmi
License: MIT
Group: SmartMet/Engines
URL: https://github.com/fmidev/smartmet-engine-contour
Source0: %{name}.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)
BuildRequires: gcc-c++
BuildRequires: make
BuildRequires: boost-devel
BuildRequires: geos-devel
BuildRequires: gdal-devel
BuildRequires: libconfig-devel
BuildRequires: smartmet-library-spine-devel >= 17.9.8
BuildRequires: smartmet-library-newbase-devel >= 17.9.11
BuildRequires: smartmet-library-macgyver-devel >= 17.8.28
BuildRequires: smartmet-library-tron >= 17.9.12
BuildRequires: smartmet-library-gis-devel >= 17.9.12
BuildRequires: sparsehash-devel
Requires: smartmet-library-gis >= 17.9.12
Requires: geos
Requires: gdal
Requires: libconfig
Requires: smartmet-library-newbase >= 17.9.11
Requires: smartmet-library-macgyver >= 17.8.28
Requires: smartmet-library-spine >= 17.9.8
%if 0%{rhel} >= 7
Requires: boost-date-time
Requires: boost-filesystem
Requires: boost-iostreams
Requires: boost-regex
Requires: boost-system
Requires: boost-thread
%endif
Provides: %{SPECNAME}
Obsoletes: smartmet-brainstorm-contour < 16.11.1
Obsoletes: smartmet-brainstorm-contour-debuginfo < 16.11.1

%description
SmartMet contour engine

%package -n %{SPECNAME}-devel
Summary: SmartMet %{SPECNAME} development headers
Group: SmartMet/Development
Provides: %{SPECNAME}-devel
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
