PROG = $(patsubst %.cpp,%,$(wildcard *Test.cpp))

REQUIRES = geos gdal configpp

include $(shell echo $${PREFIX-/usr})/share/smartmet/devel/makefile.inc

CFLAGS = -DUNIX -DUSE_UNSTABLE_GEOS_CPP_API -O0 -g $(FLAGS)

INCLUDES += \
	-I../contour \
	-isystem /usr/include/smartmet

LIBS += \
	../contour.so \
	/usr/share/smartmet/engines/querydata.so \
	$(REQUIRED_LIBS) \
	-lsmartmet-timeseries \
	-lsmartmet-spine \
	-lsmartmet-macgyver \
	-lsmartmet-newbase \
	-lsmartmet-gis \
	-lsmartmet-trax \
	-lboost_date_time \
	-lboost_thread \
	-lboost_filesystem \
	-lboost_regex \
	-lboost_iostreams \
	-lboost_timer \
	-lboost_system \
	-lbz2 -lz -lpthread

all: $(PROG)
clean:
	rm -f $(PROG) *~

test: $(PROG)
	@echo Running tests:
	@for prog in $(PROG); do ./$$prog; done

$(PROG) : % : %.cpp ../contour.so
	$(CXX) $(CFLAGS) -o $@ $@.cpp $(INCLUDES) $(LIBS)
