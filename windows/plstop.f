c$Id: plstop.f,v 1.1 2006/11/20 20:34:56 rlt Exp $
      subroutine plstop()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Close any open plot windows and stop execution

c      Inputs:
c         none

c      Outputs:
c         none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'allotd.h'
      include  'allotn.h'
      include  'counts.h'
      include  'endata.h'
      include  'iofile.h'
      include  'pdatps.h'
      include  'plflag.h'
      include  'rdata.h'
      include  'rdat1.h'
      include  'tdata.h'

      logical   lopen, setval,palloc
      character dd*5
      integer   status,vclrwk,vclswk, m,mm, ip,pr
      real*4    etime,tary(2)

      save

c     Close PostScript file if open

      if (hdcpy) call fpplcl()

      if(everon) then
        status = vclrwk()
        status = vclswk()
      endif

c     Write last log record for transient problems

      inquire(unit = ilg, opened=lopen)
      if(lopen) then
        if(niter.gt.1) then
          write(ilg,3001) nstep,niter,nform,ttim,dt,rel0,rnorm,rnmax,
     &                    aengy,etime(tary)
        else
          write(ilg,3002) nstep,niter,nform,ttim,dt,rel0,rnmax,
     &                    etime(tary)
        endif
        titer = titer + niter
        tform = tform + nform
        write(ilg,3003) titer,tform
        close(unit=ilg)
      endif

      rfl = .false.
      call ptimpl()

c     Delete existing files

      call pdelfl()

c     Close open tplot files (keep)

      do m = 30,75
        inquire(unit = m, opened = lopen)
        if(lopen) then
          close(m)
        endif
      end do ! m

c     Close last output file

      close(unit=iow)

c     Set screen back to original setting

      call pscreen()

c     Delete memory use

      mm = ndict
      do m = mm,1,-1
        ip     = dlist(m)
        dd     = dict(m)
        pr     = iprec(m)
        setval = palloc(ip,dd,0,pr)
      end do ! m

c     Stop all execution

      stop

c     Format

3001  format(2i6,i5,1p,1e11.3,1p,5e10.2,    0p,1f10.2)
3002  format(2i6,i5,1p,1e11.3,1p,2e10.2,10x,1p,1e10.2,
     &       10x,0p,1f10.2)
3003  format(/' Total',i6,i5/'  ',34('-'),' END OF FEAP LOG ',35('-'))

      end
