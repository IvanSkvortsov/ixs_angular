CXX=g++ -std=c++11
CXXFLAGS=-I .. -I ../mmap -I ../math -I ../memorystream -I ../matrix -I ../alpha 

memory_obj=../memorystream/memorystream.o ../mmap/memory.map.o
alpha_map_obj=../alpha/alpha.map.o

libraries=-L../mmap -lmmap -L../math/ -lmath -lgmp -lmpfr
.PHONY: clean

ixs.angular.map.o: ixs.angular.map.h ixs.angular.map.hpp ixs.angular.map.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ ixs.angular.map.cpp

ixs.angular.map.exe: ixs.angular.map.demo.h ixs.angular.map.main.cpp ixs.angular.map.o
	$(CXX) $(CXXFLAGS) -o $@ ixs.angular.map.main.cpp ixs.angular.map.o $(alpha_map_obj) $(memory_obj) $(libraries)

clean:
	rm *.o
