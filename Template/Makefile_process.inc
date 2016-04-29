LIBDIR=../../lib
INCDIR=../../src
MODELLIB=model_%(model)s
CXXFLAGS= -std=c++11 -O3 -fcx-fortran-rules -fcx-limited-range
LIBFLAGS= -L$(LIBDIR) -l$(MODELLIB)

all: %(process_class)s.o

%(process_class)s.o: %(process_class)s.cc $(LIBDIR)/lib$(MODELLIB).a
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(LIBFLAGS) -I$(INCDIR)

$(LIBDIR)/lib$(MODELLIB).a:
	pushd ../../src && make && popd

.PHONY: clean

clean:
	rm -f *.o
	pushd ../../src && make clean && popd
