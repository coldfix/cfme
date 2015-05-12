CFLAGS=-std=c++11 -O3
LFLAGS=-lglpk -lboost_system -lboost_timer


all: bin/redundancy bin/init-cca bin/eliminate bin/next-layer


bin/%: %.o fm.o util.o
	./generate_git_info.sh >git_info.cpp
	g++ $(LFLAGS) $^ git_info.cpp -o $@

%.o: %.cpp
	g++ $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o bin/*
