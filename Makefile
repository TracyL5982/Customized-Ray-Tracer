all: main.cpp run

SOURCES    = main.cpp
OBJECTS    = $(SOURCES:.cpp=.o)

.cpp.o:
	g++ -w -O3 -c $< -o $@

run: main.o
	g++ $(OBJECTS) $(LDFLAGS) -o $@

clean:
	rm -f *.o

