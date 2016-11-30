PROG = $(patsubst %.cpp,%,$(wildcard *Test.cpp))

MAINFLAGS = -std=c++0x -Wall -W -Wno-unused-parameter

EXTRAFLAGS = -Werror -pedantic -Wpointer-arith -Wcast-qual \
	-Wcast-align -Wconversion -Winline \
	-Wctor-dtor-privacy -Wnon-virtual-dtor -Wno-pmf-conversions \
	-Wsign-promo -Wchar-subscripts -Wold-style-cast \
	-Wshadow -Woverloaded-virtual

DIFFICULTFLAGS = -Weffc++ -Wunreachable-code -Wredundant-decls

CC = g++
CFLAGS = -DUNIX -O0 -g $(MAINFLAGS) -Wno-write-strings

INCLUDES = \
	-I ../include \
	-I /usr/include/smartmet \
	-I /usr/include/smartmet/newbase \
	`pkg-config --cflags libconfig++`

LIBS = \
	../contour.so \
	/usr/share/smartmet/engines/querydata.so \
	-lsmartmet-spine \
	-lsmartmet_gis \
	-lsmartmet_newbase \
	-lgeos \
	-lgdal \
	-lboost_date_time \
	-lboost_iostreams \
	-lboost_filesystem \
	-lboost_regex \
	-lboost_thread \
	-lboost_timer \
	-lboost_system \
	`pkg-config --libs libconfig++` \
	-lbz2 -ljpeg -lpng -lz

all: $(PROG)
clean:
	rm -f $(PROG) *~

test: $(PROG)
	@echo Running tests:
	@for prog in $(PROG); do \
	./$$prog; \
	done

$(PROG) : % : %.cpp ../contour.so
	$(CC) $(CFLAGS) -o $@ $@.cpp $(INCLUDES) $(LIBS)
