CXX = g++
FLAGS = -O2
INC = -I./primesieve/include/
LIBS = -lgmp -lgmpxx
LIBS2 = -lpthread

PRIMESIEVE_OBJS = ./primesieve/src/Erat.o ./primesieve/src/EratBig.o ./primesieve/src/EratMedium.o ./primesieve/src/EratSmall.o ./primesieve/src/PreSieve.o \
   ./primesieve/src/CpuInfo.o ./primesieve/src/MemoryPool.o ./primesieve/src/PrimeGenerator.o ./primesieve/src/PrimeSieve.o \
   ./primesieve/src/IteratorHelper.o ./primesieve/src/LookupTables.o ./primesieve/src/popcount.o ./primesieve/src/nthPrime.o ./primesieve/src/PrintPrimes.o \
   ./primesieve/src/ParallelSieve.o ./primesieve/src/iterator.o ./primesieve/src/api.o ./primesieve/src/SievingPrimes.o

.PHONY: all

all: powerAbundance powerTrialFactoring verifyPrimePowerAbundance

powerAbundance: powerAbundance.o
	$(CXX) -o $@ $^ $(LIBS)

powerTrialFactoring: powerTrialFactoring.o $(PRIMESIEVE_OBJS)
	$(CXX) -o $@ $^ $(LIBS) $(LIBS2)

verifyPrimePowerAbundance: verifyPrimePowerAbundance.o
	$(CXX) -o $@ $^ $(LIBS)

%.o: %.cpp
	$(CXX) $(FLAGS) $(INC) -c -o $@ $<

clean:
	rm -f powerAbundance powerTrialFactoring verifyPrimePowerAbundance *.o ./primesieve/src/*.o
