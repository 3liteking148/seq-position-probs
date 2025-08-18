CXXFLAGS = -Wall -O3 -g

all: bin/dummer bin/dummerl

bin/dummer: dummer.cc dummer-util.hh version.hh
	mkdir -p bin
	$(CXX) -DDOUBLE $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ dummer.cc

bin/dummerl: dummer.cc dummer-util.hh version.hh
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ dummer.cc

clean:
	rm -f bin/dummer*

VERSION1 = git describe --dirty
VERSION2 = echo '$Format:%d$ ' | sed -e 's/.*tag: *//' -e 's/[,) ].*//'

VERSION = \"`test -e .git && $(VERSION1) || $(VERSION2)`\"

version.hh: FORCE
	echo $(VERSION) | cmp -s $@ - || echo $(VERSION) > $@

FORCE:
