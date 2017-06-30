
#PROGNAME = bremsstrahlung
PROGNAME = bethecross
#PROGNAME = tsai_born_int

EXECNAME = Output

CPP = g++ -g -O2 -Wall
#CPP = g++


#INCLUDES = -I $(ROOTSYS)/include 
# (don't need if it's in CPLUS_INCLUDE_PATH)

LIBS = -lgsl -lblas
ROOTLIBS = `root-config --cflags --libs`

$(EXECNAME): $(PROGNAME).o 
	$(CPP) $(ROOTLIBS) $(LIBS) -o $@ $<

.cpp.o:
	$(CPP) $(INCLUDES) -c -o $@ $<

