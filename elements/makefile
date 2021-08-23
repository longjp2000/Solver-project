archive: 
	(cd frame; make Feap)
	(cd material; make Feap)
	(cd material/fmall; make Feap)
	(cd material/finite; make Feap)
	(cd shells; make Feap)
	(cd solid1d; make Feap)
	(cd solid2d; make Feap)
	(cd solid3d; make Feap)
	(cd thermal; make Feap)

clean: 
	(cd frame; make clean)
	(cd material; make clean)
	(cd material/small; make clean)
	(cd material/finite; make clean)
	(cd shells; make clean)
	(cd solid1d; make clean)
	(cd solid2d; make clean)
	(cd solid3d; make clean)
	(cd thermal; make clean)
	rcsclean

checkout:
	co -q frame/makefile
	(cd frame; make checkout)
	co -q material/makefile
	(cd material; make checkout)
	co -q material/small/makefile
	(cd material/small; make checkout)
	co -q material/finite/makefile
	(cd material/finite; make checkout)
	co -q shells/makefile
	(cd shells; make checkout)
	co -q solid1d/makefile
	(cd solid1d; make checkout)
	co -q solid2d/makefile
	(cd solid2d; make checkout)
	co -q solid3d/makefile
	(cd solid3d; make checkout)
	co -q thermal/makefile
	(cd thermal; make checkout)

install: 
	(cd frame; make install)
	(cd material; make install)
	(cd material/small; make install)
	(cd material/finite; make install)
	(cd shells; make install)
	(cd solid1d; make install)
	(cd solid2d; make install)
	(cd solid3d; make install)
	(cd thermal; make install)

rcs:
	(cd frame; make rcs)
	(cd material; make rcs)
	(cd material/small; make rcs)
	(cd material/finite; make rcs)
	(cd shells; make rcs)
	(cd solid1d; make rcs)
	(cd solid2d; make rcs)
	(cd solid3d; make rcs)
	(cd thermal; make rcs)
	ci -t-"" makefile

