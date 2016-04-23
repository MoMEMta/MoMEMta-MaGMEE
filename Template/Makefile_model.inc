LIBDIR=../lib
CXXFLAGS= -O3 -std=c++11 -I. -fcx-fortran-rules -fcx-limited-range

target=$(LIBDIR)/libmodel_%(model)s.a

all: $(target)

objects=HelAmps_%(model)s.o Parameters_%(model)s.o read_slha.o

$(target): $(objects)
	$(AR) cru $@ $(objects)
	ranlib $(target)

.PHONY: clean

clean:
	rm -f $(target)
	rm -f $(objects)

