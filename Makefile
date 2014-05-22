CXX = g++
CFLAGS = -Wall -O3
LDFLAGS = -lgsl -lcblas
LIBPATH = -L/opt/local/lib
LIB = Matrix.o util.o ReadSimulator.o
HEADERS = *.h

%.o:	%.cpp $(HEADERS)
	$(CXX) -c $(CFLAGS) $< -o $*.o

all: $(LIB)  SeqAssembly FastqSimulate

clean:
	rm -f *.o *~ SeqAssembly FastqSimulate

SeqAssembly: $(LIB) SeqAssembly.o
	$(CXX) $(CFLAGS) $(LIBPATH) -o SeqAssembly SeqAssembly.o $(LIB) $(LDFLAGS)

FastqSimulate: $(LIB) FastqSimulate.o
	$(CXX) $(CFLAGS) $(LIBPATH) -o FastqSimulate FastqSimulate.o $(LIB) $(LDFLAGS)
