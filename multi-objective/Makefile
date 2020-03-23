all: search

search: main.o test_problem.o
	g++ -ggdb -o search main.o test_problem.o

main.o: main.cpp test_problem.cpp
	g++ -ggdb -c main.cpp test_problem.cpp

dep:
	echo "Do nothing"

clean:
	rm -f main *.o

