c$Id: ploadc.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine ploadc()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Control routine for data inputs based on nodal
c               coordinates

c      Inputs:
c         none      - Data retrieved through common blocks

c      Outputs:
c         none      - Data stored in blank common locations
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'corset.h'
      include  'corfil.h'
      include  'print.h'
      include  'sdata.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   setvar,palloc
      integer   ityp(11), ndft(11), n,nn

      save

c     Allocate temporary arrays

      setvar = palloc(111,'TEMP1',numnp+4*numel,1)
      setvar = palloc(112,'TEMP2',numnp        ,2)

c     Set the counter

      nn = 0

c     Set coordinate angle conditions

      if(angfl) then
        nn       = nn + 1
        ndft(nn) = 1
        ityp(nn) = 3
      endif

c     Set coordinate boundary codes

      if(boufl) then
        nn       = nn + 1
        ndft(nn) = ndf
        ityp(nn) = 2
      endif

c     Set cordinate force values

      if(forfl) then
        nn       = nn + 1
        ndft(nn) = ndf
        ityp(nn) = 5
      endif

c     Set cordinate displacement values

      if(disfl) then
        nn       = nn + 1
        ndft(nn) = ndf
        ityp(nn) = 4
      endif

c     Set cordinate force proportional load values

      if(cprfl) then
        nn       = nn + 1
        ndft(nn) = ndf
        ityp(nn) = 6
      endif

c     Set coordinate surface loads

      if(surfl) then
        nn       = nn + 1
        ndft(nn) = ndf
        ityp(nn) = 1
      endif

c     Set cordinate force proportional load values

      if(damfl) then
        nn       = nn + 1
        ndft(nn) = ndf
        ityp(nn) = 7
      endif

c     Set cordinate force proportional load values

      if(masfl) then
        nn       = nn + 1
        ndft(nn) = ndf
        ityp(nn) = 8
      endif

c     Set cordinate force proportional load values

      if(stifl) then
        nn       = nn + 1
        ndft(nn) = ndf
        ityp(nn) = 9
      endif

c     Set coordinate boundary codes

      if(basfl) then
        nn       = nn + 1
        ndft(nn) = ndf
        ityp(nn) = 10
      endif

c     Set cordinate force proportional load values

      if(lfrfl) then
        nn       = nn + 1
        ndft(nn) = ndf
        ityp(nn) = 11
      endif

c     Call routine to do generations

      do n = 1,nn
        call pesurf(mr(np(31)),mr(np(33)),mr(np(111)),mr(np(111)+numnp),
     &              hr(np(112)),hr(np(43)),hr(np(27)),hr(np(45)),
     &              ndft(n),ndm,nen,nen1,numnp,numel,prt,prth,ityp(n))
      end do ! n

c     Destroy temporary arrays

      setvar = palloc(112,'TEMP2',0,2)
      setvar = palloc(111,'TEMP1',0,1)

      end
