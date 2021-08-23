checkout:   
	co -q RCS/*.f,v
	@@echo "--> HP checked out <--"

clean:
	rcsclean -q
	rm -f *.o
	@@echo "--> HP rcscleaned <--"

rcs:
	mkdir RCS
	ci -t-"" *.f
	ci -t-"" makefile
	@@echo "--> HP created RCS<--"

