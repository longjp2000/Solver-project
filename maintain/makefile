include ../makefile.in

clean:
	rcsclean -q
	rm -f *.o

checkout:
	co -q RCS/*.c,v

rcs:
	mkdir RCS
	ci -t-"" *.c
	ci -t-"" makefile

