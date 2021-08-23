c$Id: zzprod.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine zzprod(lct,ctl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Driver routine for Zienkiewicz-Zhu stress projection

c      Inputs:
c         lct     - Character array for options
c                   'off' - return to lumped nodal projection
c         ctl     - Real parameter for material number to project

c      Outputs:
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      logical    pcomp,setvar,palloc
      character  lct*15
      integer    ma,numst,ipmax
      real*8     ctl

      include   'bdata.h'
      include   'cdata.h'
      include   'fdata.h'
      include   'iofile.h'
      include   'pbody.h'
      include   'pdata3.h'
      include   'prstrs.h'
      include   'sdata.h'

      include   'pointer.h'
      include   'comblk.h'

      save

c     Reset to lumped projection

      if(pcomp(lct,'off ',4)) then

        fl(11) = .false.
        if(ior.lt.0) then
          write(*,2000)
        endif
        write(iow,2000)

c     Assemble patches for ZZ-projection

      else

c       Set parameters

        numst = npstr - 1
        ma    = nint(ctl)
        maplt = ma
        if(ma.eq.0) then
          if(ior.lt.0) then
            write(*,2001)
          endif
          write(iow,2001)
        else
          if(ior.lt.0) then
            write(*,2002) ma
          endif
          write(iow,2002) ma
        endif

c       Set arrays to find patches

        if (plfl) then
          setvar = palloc( 58,'NDNP ',numnp*npstr, 2)
          setvar = palloc( 57,'NDER ',numnp*8    , 2)
          setvar = palloc( 60,'NDNS ',max(nen,npstr,nst,nst),2)
          setvar = palloc(207,'NSCR ',numel      , 2)
          plfl = .false.
        endif
        nph  = np( 58)
        ner  = np( 57)

        setvar = palloc(218,'ZZIB ',numnp    , 1)
        setvar = palloc(219,'ZZIP ',numnp    , 1)

c       Call routines to compute Zienkiewicz-Zhu projection

        call zzpro1(mr(np(33)),mr(np(218)),mr(np(219)),
     &              ma,ndm,nen,nen1,numnp,numel,1,ipmax)

        setvar = palloc(111,'TEMP1 ',ipmax, 1)

        call zzpro2(mr(np(33)),mr(np(218)),mr(np(219)),mr(np(111)),
     &              hr(np(43)),hr(nph+numnp),hr(ner+numnp),
     &              ndm,nen,nen1,numnp,numel,numst)

        setvar = palloc(111,'TEMP1 ',0, 1)

        fl(11) = .true.

      endif

c     Formats

2000  format(/'      Disable ZZ-projections - Use lumped method'/)
2001  format(/'      Perform ZZ-projections for all materials'/)
2002  format(/'      Perform ZZ-projections for material',i4/)

      end
