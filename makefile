compile: proj.cpp
	g++ -c proj.cpp
	g++ -o proj proj.cpp
	./proj
	rm *.o
