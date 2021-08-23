# The place where you need to make selections is under target
# "install".  Under this target only, appropriately comment or
# uncomment the line associated with your Fortran compiler type
# (i.e., f77 or f90).

# N.B.  It is necessary to also modify 'makefile.in' appropriately.
feap:
	(cd main; make feap)
	@@echo "--> Feap executable made <--"

#archive:   
#	(cd program; make archive)
#	(cd contact; make archive)
#	(cd elements; make archive)
#	(cd plot; make archive)
#	(cd unix; make archive)
#	(cd user; make archive)
#	@@echo "--> Feap Archive updated <--"
#
#checkout:   
#	co -q makefile
#	co -q makefile.in
#	co -q maintain/makefile
#	(cd maintain; make checkout)
#	co -q program/makefile
#	(cd program; make checkout)
#	co -q contact/makefile
#	(cd contact; make checkout)
#	co -q elements/makefile
#	(cd elements; make checkout)
#	co -q memory/makefile
#	(cd memory; make checkout)
#	co -q plot/makefile
#	(cd plot; make checkout)
#	co -q unix/makefile
#	(cd unix; make checkout)
#	co -q user/makefile
#	(cd user; make checkout)
#	co -q main/makefile
#	(cd main; make checkout)
#	co -q windows/makefile
#	(cd windows; make checkout)
#	co -q window2/makefile
#	(cd window2; make checkout)
#	co -q hp/makefile
#	(cd hp; make checkout)
#	co -q ibm/makefile
#	(cd ibm; make checkout)
#	co -q f77/makefile
#	(cd f77; make checkout)
#	co -q f90/makefile
#	(cd f90; make checkout)
#	co -q packages/arpack/makefile
#	(cd packages/arpack; make checkout)
#	co -q packages/arpack/archive/makefile
#	(cd packages/arpack/archive; make checkout)
#	co -q packages/blas/makefile
#	(cd packages/blas; make checkout)
#	co -q packages/lapack/makefile
#	(cd packages/lapack; make checkout)
#	co -q packages/meshmod/makefile
#	(cd packages/meshmod; make checkout)
#	co -q parfeap/makefile
#	(cd parfeap; make checkout)
#	co -q parfeap/packages/arpack/makefile
#	(cd parfeap/packages/arpack; make checkout)
#	co -q parfeap/partition/makefile
#	(cd parfeap/partition;  make checkout)
#	co -q parfeap/unix/makefile
#	(cd parfeap/unix;  make checkout)
#	co -q parfeap/windows/makefile
#	(cd parfeap/windows;  make checkout)
#	co -q patch/makefile
#	(cd patch; make checkout)
#	(cd include; make checkout)
#	(cd include/integer4; make checkout)  # For integer*4 pointers
#	(cd include/integer8; make checkout)  # For integer*8 pointers
#	@@echo "--> Feap checked out <--"

clean:
	rm -f *.a
	(cd include; make clean)
	(cd include/integer4; make clean)     # For integer*4 pointers
	(cd include/integer8; make clean)     # For integer*8 pointers
	(cd maintain; make clean)
	(cd program; make clean)
	(cd contact; make clean)
	(cd elements; make clean)
	(cd memory; make clean)
	(cd plot; make clean)
	(cd unix; make clean)
	(cd user; make clean)
	(cd main; make clean)
	(cd windows; make clean)
	(cd window2; make clean)
	(cd hp; make clean)
	(cd ibm; make clean)
	(cd f77; make clean)
	(cd f90; make clean)
	#(cd packages/arpack; make clean)
	#(cd packages/arpack/archive; make clean)
	#(cd packages/blas; make clean)
	#(cd packages/lapack; make clean)
	#(cd packages/meshmod; make clean)
	#(cd parfeap; make clean)
	#(cd parfeap/packages/arpack; make clean)
	#(cd parfeap/partition;  make clean)
	#(cd parfeap/unix;  make clean)
	#(cd parfeap/windows;  make clean)
	#(cd patch; make clean)
	rcsclean -q
	@@echo "--> Feap rcscleaned <--"

install:
	(cd program; make install)
	(cd contact; make install)
	(cd elements; make install)
	(cd plot; make install)
	(cd unix; make install)
	(cd user; make install)
	(cd f77; make install)
#	(cd f90; make install)
	(cd main; make feap)
	@@echo "--> Feap Installed <--"

#rcs:
#	(mkdir RCS; ci -t-"" makefile; ci -t-"" makefile.in)
#	(cd maintain; make rcs)
#	(cd program; make rcs)
#	(cd contact; make rcs)
#	(cd elements; make rcs)
#	(cd plot; make rcs)
#	(cd unix; make rcs)
#	(cd user; make rcs)
#	(cd main; make rcs)
#	(cd memory; make rcs)
#	(cd hp; make rcs)
#	(cd ibm; make rcs)
#	(cd f77; make rcs)
#	(cd f90; make rcs)
#	(cd packages/arpack; make rcs)
#	(cd packages/blas; make rcs)
#	(cd packages/lapack; make rcs)
#	(cd packages/meshmod; make rcs)
#	(cd parfeap; make rcs)
#	(cd parfeap/packages/arpack; make rcs)
#	(cd parfeap/partition;  make rcs)
#	(cd parfeap/unix;  make rcs)
#	(cd parfeap/windows;  make rcs)
#	(cd patch; make rcs)
#	(cd windows; make rcs)
#	(cd window2; make rcs)
#	(cd include; make rcs)
#	(cd include/integer4; make rcs)
#	(cd include/integer8; make rcs)
#	@@echo "--> Feap created RCS<--"

