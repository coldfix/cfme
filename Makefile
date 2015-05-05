CFLAGS=-I$(FM)/include -std=c++11 -O3
LFLAGS=-L$(FM)/lib -lglpk


all: bin/redundancy bin/cfme


bin/redundancy: redundancy.o fm.o
	g++ $(LFLAGS) $^ -o $@

bin/cfme: cfme.o fm.o
	g++ $(LFLAGS) $^ -o $@

%.o: %.cpp
	g++ $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o bin/*
