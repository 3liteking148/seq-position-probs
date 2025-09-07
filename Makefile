CXXFLAGS = -Wall -O3 -g -march=native

all: bin/dummer bin/dummerl bin/dummer-build

bin/dummer: dummer.cc dummer-util.hh version.hh can_i_haz_simd.hh
	mkdir -p bin
	$(CXX) -DDOUBLE $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ dummer.cc

bin/dummerl: dummer.cc dummer-util.hh version.hh can_i_haz_simd.hh
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ dummer.cc

bin/dummer-build: dummer-build.cc dummer-util.hh version.hh priors/gap-priors.hh priors/wheeler4.hh priors/blocks9.hh
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ dummer-build.cc

clean:
	rm -f bin/dummer*

VERSION1 = git describe --dirty
VERSION2 = echo '$Format:%d$ ' | sed -e 's/.*tag: *//' -e 's/[,) ].*//'

VERSION = \"`test -e .git && $(VERSION1) || $(VERSION2)`\"

version.hh: FORCE
	echo $(VERSION) | cmp -s $@ - || echo $(VERSION) > $@

FORCE:
