FM=/opt/fm-lib

CFLAGS=-I$(FM)/include -std=c++11
LFLAGS=-L$(FM)/lib -lfm -Wl,-rpath=$(FM)/lib

cfme: main.o cfme.o fm.o
	g++ $(LFLAGS) $^ -o $@

%.o: %.cpp
	g++ $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o
