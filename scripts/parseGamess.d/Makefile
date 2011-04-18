all: parseGamess

parseGamess: parseGamess.o
	g++ -o parseGamess parseGamess.o
	cp parseGamess ../

parseGamess.o: parseGamess.cpp
	g++ -o parseGamess.o -c parseGamess.cpp

clean:
	rm parseGamess *.o
