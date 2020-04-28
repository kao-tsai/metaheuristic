all: search

search: main.o test_problem.o
	g++ -ggdb -O3 -o search main.o test_problem.o

main.o: main.cpp
	g++ -ggdb -O3 -c main.cpp

test_problem.o: test_problem.cpp
	g++ -ggdb -O3 -c test_problem.cpp

dep:
	echo "Do nothing"

clean:
	rm -f main *.o

