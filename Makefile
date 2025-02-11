CXXFLAGS = -Wall -O3

bin/seq-position-probs: seq-position-probs.cc
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ seq-position-probs.cc
