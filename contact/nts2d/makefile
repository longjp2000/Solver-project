include ../../makefile.in

SOURCES = $(FSOURCE)*.$(FEXT)

archive: $(SOURCES)
	$(AR) $(ARFEAP) $%
	@@echo "--->Contact NTS2D objects now up-to-date<---"
	rm -f $% *.o

clean:
	rcsclean -q
	rm -f *.o

checkout:
	co -q RCS/*.f,v

install: $(SOURCES)
	$(FF) -c -I$(FINCLUDE) $(FFOPTFLAG) *.f
	$(AR) $(ARFEAP) *.o
	@@echo "--->Contact NTS2D objects now installed in archive<---"
	rm -f $% *.o

rcs:
	mkdir RCS
	ci -t-"" *.f
	ci -t-"" makefile

.f.o:
	$(FF) -c -I$(FINCLUDE) $(FFOPTFLAG) *.f -o $*.o

.c.o:
	$(CC) -c $(CCOPTFLAG) *.c -o $*.o

