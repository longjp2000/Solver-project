checkout:   
	co -q RCS/*.f,v
	@@echo "--> IBM checked out <--"

clean:
	rcsclean -q
	rm -f *.o
	@@echo "--> IBM rcscleaned <--"

rcs:
	mkdir RCS
	ci -t-"" *.f
	ci -t-"" makefile
	@@echo "--> IBM created RCS<--"

