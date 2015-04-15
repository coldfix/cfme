FM=/opt/fm-lib

CFLAGS=-I$(FM)/include -std=c++11
LFLAGS=-L$(FM)/lib -lfm -lglpk -Wl,-rpath=$(FM)/lib


all: redundancy cfme


redundancy: redundancy.o fm.o
	g++ $(LFLAGS) $^ -o $@

cfme: cfme.o fm.o
	g++ $(LFLAGS) $^ -o $@

%.o: %.cpp
	g++ $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o
