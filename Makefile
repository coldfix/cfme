CFLAGS=-I$(FM)/include -std=c++11 -O3
LFLAGS=-L$(FM)/lib -lglpk


all: redundancy cfme


redundancy: redundancy.o fm.o
	g++ $(LFLAGS) $^ -o $@

cfme: cfme.o fm.o
	g++ $(LFLAGS) $^ -o $@

%.o: %.cpp
	g++ $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o
