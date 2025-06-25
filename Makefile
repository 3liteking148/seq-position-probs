CXXFLAGS = -Wall -O3 -g

bin/seq-position-probs: seq-position-probs.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ seq-position-probs.cc

clean:
	rm -f bin/seq-position-probs
