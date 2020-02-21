all: search

search: main.o
	g++ -ggdb -O3 -o search main.o 

main.o: main.cpp
	g++ -ggdb -O3 -c main.cpp

dep:
	echo "Do nothing"

clean:
	rm -f main *.o
