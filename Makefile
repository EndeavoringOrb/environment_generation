CXX=g++
CXXFLAGS=-IC:/Users/aaron/CODING/C++_INCLUDES/SFML-2.5.1/include -DSFML_STATIC
LDFLAGS=-LC:/Users/aaron/CODING/C++_INCLUDES/SFML-2.5.1/lib
LDLIBS=-lsfml-graphics-s -lsfml-window-s -lsfml-system-s -lopengl32 -lfreetype -lwinmm -lgdi32 -lsfml-main

all: make_blocks

make_blocks: make_blocks.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

make_blocks.o: make_blocks.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f make_blocks.o make_blocks