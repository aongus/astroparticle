
PROGNAME = energyintegral

EXECNAME = Output

CPP = g++ -g -O2 -Wall

#INCLUDES = -I $(ROOTSYS)/include 
# (don't need if it's in CPLUS_INCLUDE_PATH)

LIBS = -lgsl
ROOTLIBS = `root-config --cflags --libs`

$(EXECNAME): $(PROGNAME).o
	$(CPP) $(ROOTLIBS) $(LIBS) -o $@ $^

.cpp.o:
	$(CPP) $(INCLUDES) -c -o $@ $<

