include ../makefile.in

OBJECTS = feap81.o

SOURCES = $(FSOURCE)feap81.$(FEXT)

feap: $(OBJECTS) $(SOURCES) $(ARFEAP)
	ranlib $(ARFEAP)
	$(FF) -o feap $(OBJECTS) $(ARFEAP) $(LDOPTIONS)

.f.o:
	$(FF) -c $(FFOPTFLAG) -I$(FINCLUDE) $*.f -o $*.o

.c.o:
	$(CC) -c $(CCOPTFLAG)  $*.c -o $*.o

clean:
	rcsclean -q
	rm -f *.o

checkout:
	co -q RCS/*.f,v

rcs:
	mkdir RCS
	ci -t-"" *.f
	ci -t-"" makefile

