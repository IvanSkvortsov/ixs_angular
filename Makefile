CXX=g++ -std=c++11
CXXFLAGS=-I .. -I ../mmap -I ../math -I ../memorystream -I ../matrix -I ../alpha -I ../omega -I ../vector_3d_ -I ../poly_n -I ../array \
	 -I ../xyz_pown_integral -I ../omega

memory_obj=../memorystream/memorystream.o ../mmap/memory.map.o
alpha_map_obj=../alpha/alpha.map.o
matrix_slm_obj=../omega/matrix.slm.o ../poly_n/poly.n.o ../array/array.n.o
vector_3d_obj=../vector_3d_/vector.3d.o
objects=geom.slm.o ixs.angular.map.o ixs.angular.dat.o

libraries=-L../mmap -lmmap -L../math/ -lmath -lgmp -lmpfr
.PHONY: clean

ixs.angular.dat.o: ixs.angular.dat.h ixs.angular.dat.hpp ixs.angular.dat.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ ixs.angular.dat.cpp

ixs.angular.dat.exe: ixs.angular.dat.demo.h ixs.angular.dat.main.cpp $(objects) ixs.angular.mem.h
	$(CXX) $(CXXFLAGS) -o $@ ixs.angular.dat.main.cpp $(objects) $(alpha_map_obj) $(memory_obj) $(libraries)

ixs.angular.map.o: ixs.angular.map.h ixs.angular.map.hpp ixs.angular.map.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ ixs.angular.map.cpp

ixs.angular.map.exe: ixs.angular.map.demo.h ixs.angular.map.main.cpp ixs.angular.map.o ixs.angular.mem.h
	$(CXX) $(CXXFLAGS) -o $@ ixs.angular.map.main.cpp ixs.angular.map.o $(alpha_map_obj) $(memory_obj) $(libraries)

geom.slm.o: geom.slm.h geom.slm.hpp geom.slm.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ geom.slm.cpp

geom.slm.exe: geom.slm.o geom.slm.demo.h geom.slm.main.cpp
	$(CXX) $(CXXFLAGS) -o $@ geom.slm.main.cpp geom.slm.o $(vector_3d_obj) $(matrix_slm_obj) $(libraries)

clean:
	rm *.o
