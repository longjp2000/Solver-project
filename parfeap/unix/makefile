include ../../makefile.in

SOURCES = $(FSOURCE)*.f,v

OBJECTS = $(FSOURCE)*.o

feap: $(SOURCES) $(OBJECTS)
	-${FLINKER} ${FFLAGS} $*.f

clean:
	rcsclean -q
	rm -f *.o

checkout:
	co -q RCS/*.f,v

rcs:
	mkdir RCS
	ci -t-"" *.f
	ci -t-"" makefile

.f.o:
	$(FF) -c -I$(FINCLUDE) $(FFOPTFLAG) *.f -o $*.o

.c.o:
	$(CC) -c $(CCOPTFLAG) *.c -o $*.o
