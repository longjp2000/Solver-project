c$Id: flaechen.f,v 1.1 2006/11/20 20:32:46 rlt Exp $
      subroutine flaechen(cs0,ics,surpoin)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Anna Haraldsson             February 1998            1.0

c      Acronym:

c      Purpose: Get patch-elements of each node of all surfaces of "ics"
c               and store them in knotn.

c      Inputs:
c         cs0(*)  - Surface parameters
c         ics(*)  - Surface facets
c         surpoin - Surface points

c      Outputs:
c         knotn(*)- Facets connected to nodes through pointer np(191)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'c_0.h'
      include   'c_contac.h'
      include   'c_comnd.h'
      include   'c_dict.h'
      include   'pointer.h'
      include   'comblk.h'

      logical    setvar,palloc
      integer    surf,ofsurf,dnope,neps,nope, i, lockn,laenge
      integer    ics(*),surpoin(*)
      real*8     cs0(nr0,n0c1:nc01,*)

      save

      laenge = 0

      do surf = 1,numcs
        ofsurf        = cs0(2,-1,surf)
        nope          = cs0(2, 0,surf)
        dnope         = cs0(4,-1,surf)
        neps          = cs0(3,-1,surf)
        surpoin(surf) = laenge

c       Dynamic memory allocation and setting of knotn array

        setvar = palloc(136,'CTEM1',5*neps*dnope,1)
        lockn  = 1
        call einfl(ics(ofsurf),mr(np(136)),lockn,laenge,dnope,neps,nope)
        setvar = palloc( 191,'KNOTN',laenge+lockn, 1)
        do i = 0,lockn-1
          mr(np(191)+laenge+i) = mr(np(136)+i)
        end do ! i
        laenge = laenge + lockn
        setvar = palloc( 136,'CTEM1',   0, 1)

      end do ! surf

      end
