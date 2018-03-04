###### Definition of Macros ######


PROGRAM  = hplus.exe


CXX=`root-config --cxx`
CXXFLAGS=`root-config --cflags`
LDFLAGS=`root-config --ldflags`
LDLIBS=`root-config --glibs`
ROOTLIBS='-lRooFit -lHtml -lMinuit -lRooFitCore -lRooStats -lHistFactory'

#CXXFLAGS = $(shell root-config --cflags)
#ROOTLIBS = $(shell root-config --glibs)
LIBS     = $(ROOTLIBS)
#OPT      = -Wall -E -MM

GCC      = /usr/bin/gcc

OBJECTS  = main.o hplus.o model.o BR.o limit.o

###### Build Programm ######

all: $(OBJECTS)
	g++ -g  -Wall `root-config --cflags --libs` -L$(ROOTSYS)/lib $? -o $(PROGRAM)

###### Objects ######

main.o: main.cxx hplus.h model.h BR.h limit.h
	g++ -g  -Wall `root-config --cflags --libs` -L$(ROOTSYS)/lib -c main.cxx
hplus.o: hplus.cxx hplus.h
	g++ -g  -Wall `root-config --cflags --libs` -L$(ROOTSYS)/lib -c hplus.cxx
model.o: model.cxx model.h
	g++ -g  -Wall `root-config --cflags --libs` -L$(ROOTSYS)/lib -c model.cxx
BR.o: BR.cxx BR.h
	g++ -g  -Wall `root-config --cflags --libs` -L$(ROOTSYS)/lib -c BR.cxx
limit.o: limit.cxx limit.h
	g++ -g  -Wall `root-config --cflags --libs` -L$(ROOTSYS)/lib -c limit.cxx


###### make clean ######

clean:
	rm -f *.o
	rm -f $(PROGRAM)
	rm -rf ti_files
	rm -f *~;

