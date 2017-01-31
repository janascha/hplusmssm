###### Definition of Macros ######

PROGRAM  = hplus.exe

CXXFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --glibs)

LIBS     = $(ROOTLIBS)
#OPT      = -Wall -E -MM

OPT      =
GCC      = /usr/bin/gcc

OBJECTS  = main.o hplus.o model.o BR.o limit.o

###### Build Programm ######

all: $(OBJECTS)
	$(GCC) $(OPT) $(CXXFLAGS) -g $(OBJECTS) -o $(PROGRAM) $(LIBS) -lRooFit -lMathCore -lTreePlayer -lTMVA -lMLP

###### Objects ######

main.o: main.cxx hplus.h model.h BR.h limit.h
	$(GCC) $(OPT) $(CXXFLAGS) -c main.cxx
hplus.o: hplus.cxx hplus.h
	$(GCC) $(OPT) $(CXXFLAGS) -c hplus.cxx
model.o: model.cxx model.h
	$(GCC) $(OPT) $(CXXFLAGS) -c model.cxx
BR.o: BR.cxx BR.h
	$(GCC) $(OPT) $(CXXFLAGS) -c BR.cxx
limit.o: limit.cxx limit.h
	$(GCC) $(OPT) $(CXXFLAGS) -c limit.cxx


###### make clean ######

clean:
	rm -f *.o
	rm -f $(PROGRAM)
	rm -rf ti_files
	rm -f *~;

