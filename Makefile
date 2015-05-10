CFLAGS=-std=c++11 -O3
LFLAGS=-lglpk -lboost_system -lboost_timer


all: bin/redundancy bin/init-cca bin/eliminate


bin/%: %.o fm.o util.o
	g++ $(LFLAGS) $^ -o $@

%.o: %.cpp
	g++ $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o bin/*
