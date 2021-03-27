CXX = g++
FLAGS = -O2
LIBS = -lgmp -lgmpxx

.PHONY: all

all: powerAbundance verifyPrimePowerAbundance

powerAbundance: powerAbundance.o
	$(CXX) -o $@ $^ $(LIBS)

verifyPrimePowerAbundance: verifyPrimePowerAbundance.o
	$(CXX) -o $@ $^ $(LIBS)

%.o: %.cpp
	$(CXX) -c -o $@ $< $(FLAGS)

clean:
	rm -f powerAbundance verifyPrimePowerAbundance *.o
