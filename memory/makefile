include ../makefile.in

SOURCES = $(FSOURCE)*.$(FEXT)

clean:
	rcsclean -q
	rm -f *.o

checkout:
	co -q RCS/*.f,v

compile: $(SOURCES)
	$(FF) -c -I$(FINCLUDE) $(FFOPTFLAG) *.f
	@@echo "--->MEMORY objects now compiled<---"

rcs:
	mkdir RCS
	ci -t-"" *.f
	ci -t-"" makefile

.f.o:
	$(FF) -c -I$(FINCLUDE) $(FFOPTFLAG) *.f -o $*.o

