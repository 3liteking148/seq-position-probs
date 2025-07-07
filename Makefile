CXXFLAGS = -Wall -O3 -g

all: bin/seq-position-probs bin/seq-position-probs2

bin/seq-position-probs: seq-position-probs.cc version.hh
	$(CXX) -DDOUBLE $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ seq-position-probs.cc

bin/seq-position-probs2: seq-position-probs.cc version.hh
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ seq-position-probs.cc

clean:
	rm -f bin/seq-position-probs*

VERSION1 = git describe --dirty
VERSION2 = echo '$Format:%d$ ' | sed -e 's/.*tag: *//' -e 's/[,) ].*//'

VERSION = \"`test -e .git && $(VERSION1) || $(VERSION2)`\"

version.hh: FORCE
	echo $(VERSION) | cmp -s $@ - || echo $(VERSION) > $@

FORCE:
