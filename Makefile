CFLAGS=-std=c++11 -O3
LFLAGS=-lglpk -lboost_system -lboost_timer


BIN = \
	  check_equivalence \
	  init-cca \
	  eliminate \
	  next-layer \
	  check_shift_invariance \
	  diff-systems \
	  minimize_system \
	  random \


all: $(addprefix bin/,$(BIN))


bin/%: %.o fm.o util.o
	./generate_git_info.sh >git_info.cpp
	g++ $(LFLAGS) $^ git_info.cpp -o $@

%.o: %.cpp
	g++ $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o bin/*
