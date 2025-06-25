CXXFLAGS = -Wall -O3 -g

all: bin/seq-position-probs bin/seq-position-probs2

bin/seq-position-probs: seq-position-probs.cc
	$(CXX) -DDOUBLE $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ seq-position-probs.cc

bin/seq-position-probs2: seq-position-probs.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ seq-position-probs.cc

clean:
	rm -f bin/seq-position-probs*
