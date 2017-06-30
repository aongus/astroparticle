
#PROGNAME = nucross
PROGNAME = pfunc
#PROGNAME = ic80nuarea

EXECNAME = Output

CPP = g++ -g -O2 -Wall

LHALIB = -lLHAPDF

GSLLIB = -lgsl
ROOTLIBS = `root-config --cflags --libs`

$(EXECNAME): $(PROGNAME).o 
	$(CPP) $(ROOTLIBS) $(LHALIB) $(GSLLIB) -o $@ $<

.cpp.o:
	$(CPP) -c -o $@ $<

