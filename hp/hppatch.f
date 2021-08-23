c$Id: hppatch.f,v 1.1 2006/11/20 20:33:47 rlt Exp $

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Patch routines for HP computers
c-----[--.----+----.----+----.-----------------------------------------]

        real function etime (tary)
	integer clock,itime
	real ttime, tary(2)

	itime = clock()
	ttime = float(itime)/1000000.0
	tary(2) = ttime
	tary(1) = ttime
	etime   = ttime

	end
	subroutine fdate ()
	end
	subroutine flush ()
	end
	subroutine getlog ()
	end
	subroutine lnblnk ()
	end
