PROG = $(patsubst %.cpp,%,$(wildcard *Test.cpp))

FLAGS = -std=c++11 -Wall -W -Wno-unused-parameter

EXTRAFLAGS = -Werror -pedantic -Wpointer-arith -Wcast-qual \
	-Wcast-align -Wconversion -Winline \
	-Wctor-dtor-privacy -Wnon-virtual-dtor -Wno-pmf-conversions \
	-Wsign-promo -Wchar-subscripts -Wold-style-cast \
	-Wshadow -Woverloaded-virtual

DIFFICULTFLAGS = -Weffc++ -Wunreachable-code -Wredundant-decls

ifeq ($(TSAN), yes)
  FLAGS += -fsanitize=thread
endif
ifeq ($(ASAN), yes)
  FLAGS += -fsanitize=address -fsanitize=pointer-compare -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize-address-use-after-scope
endif


CC = g++
CFLAGS = -DUNIX -O0 -g $(FLAGS) -Wno-write-strings

# Boost 1.69

ifneq "$(wildcard /usr/include/boost169)" ""
  INCLUDES += -isystem /usr/include/boost169
  LIBS += -L/usr/lib64/boost169
endif

ifneq "$(wildcard /usr/gdal30/include)" ""
  INCLUDES += -isystem /usr/gdal30/include
  LIBS += -L$(PREFIX)/gdal30/lib
else
  INCLUDES += -isystem /usr/include/gdal
endif


INCLUDES += \
	-I../contour \
	-I/usr/include/smartmet \
	-I/usr/include/smartmet/newbase \
	`pkg-config --cflags libconfig++`

LIBS += \
	../contour.so \
	/usr/share/smartmet/engines/querydata.so \
	-lsmartmet-spine \
	-lsmartmet-gis \
	-lsmartmet-newbase \
	-lgeos \
	-lgdal \
	-lboost_date_time \
	-lboost_iostreams \
	-lboost_filesystem \
	-lboost_thread \
	-lboost_timer \
	-lboost_system \
	`pkg-config --libs libconfig++` \
	-lbz2 -ljpeg -lpng -lz \
	-lpthread

all: $(PROG)
clean:
	rm -f $(PROG) *~

test: $(PROG)
	@echo Running tests:
	@for prog in $(PROG); do \
	./$$prog; \
	done

$(PROG) : % : %.cpp ../contour.so
	$(CXX) $(CFLAGS) -o $@ $@.cpp $(INCLUDES) $(LIBS)
