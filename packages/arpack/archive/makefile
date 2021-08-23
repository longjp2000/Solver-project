
include $(FEAPHOME8_1)/makefile.in

SOURCES = *.f

OBJECTS = *.o

clean: 
	rcsclean -q
	rm -f *.o

checkout:
	co -q RCS/*.f,v
	co -q RCS/*.h,v

install: $(SOURCES)
	$(FF) -c $(FFOPTFLAG) *.f
	ar -rv arpacklib.a *.o
	ranlib arpacklib.a
	@@echo "--->ARPACK objects now installed in archive<---"
	rm -f $% *.o

rcs:
	mkdir RCS
	ci -t-"" *.f
	ci -t-"" *.h
	ci -t-"" makefile

.f.o:
	$(FF) -c $(FFOPTFLAG) *.f -o $*.o

.c.o:
	$(CC) -c $(CCOPTFLAG) *.c -o $*.o

