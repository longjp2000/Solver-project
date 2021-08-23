include ../makefile.in

SOURCES = $(FSOURCE)*.$(FEXT)

archive: $(SOURCES)
	$(AR) $(ARFEAP) $%
	@@echo "--->UNIX objects now up-to-date<---"
	rm -f $% *.o

clean:
	rcsclean -q
	rm -f *.o

checkout:
	co -q RCS/*.f,v
	co -q RCS/*.c,v

install: $(SOURCES)
	$(FF) -c -I$(FINCLUDE) $(FFOPTFLAG) *.f
	$(CC) -c -I$(CINCLUDE) $(CCOPTFLAG) *.c
	$(AR) $(ARFEAP) *.o
	@@echo "--->UNIX objects now installed in archive<---"
	rm -f $% *.o

rcs:
	mkdir RCS
	ci -t-"" *.f
	ci -t-"" *.c
	ci -t-"" makefile

.f.o:
	$(FF) -c -I$(FINCLUDE) $(FFOPTFLAG) *.f -o $*.o

.c.o:
	$(CC) -c -I$(CINCLUDE) $(CCOPTFLAG) *.c -o $*.o

