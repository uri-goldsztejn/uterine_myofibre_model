CCFLAGS = -std=c++11 -Wall -I/home/research/uri.goldsztejn/gsl/include/gsl
CLFLAG = -L/home/research/uri.goldsztejn/gsl/lib

all : stringOfCells_main.o StringOfCellsContainer.o Stimulator.o Cell.o
	g++ $(CLFLAG) -o myofibreModel stringOfCells_main.o StringOfCellsContainer.o Stimulator.o Cell.o -lgsl -lgslcblas -lm 

stringOfCells_main.o : stringOfCells_main.cpp  Stimulator.h Cell.h StringOfCellsContainer.h
	g++ $(CCFLAGS) -c stringOfCells_main.cpp

StringOfCellsContainer.o : StringOfCellsContainer.cpp StringOfCellsContainer.h
	g++ $(CCFLAGS) -c StringOfCellsContainer.cpp
	
Stimulator.o : Stimulator.cpp Stimulator.h
	g++ $(CCFLAGS) -c Stimulator.cpp
	
Cell.o : Cell.cpp Cell.h
	g++ $(CCFLAGS) -c Cell.cpp

clean : all
	rm -f stringOfCells_main.o StringOfCellsContainer.o Stimulator.o Cell.o
