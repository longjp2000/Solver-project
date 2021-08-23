c$Id: setcomp.f,v 1.1 2006/11/20 20:32:46 rlt Exp $
      subroutine setcomp (npair,cs0,cm0,cp0,hic)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: SET COMmons for a contact Pair

c      Purpose: Update values related to the contactpair in commons

c      Inputs :
c         npair   - # of contact pair
c         cs0(*)  - Contact surfaces control data
c         cm0(*)  - Contact material control data
c         cp0(*)  - Contactpair control data
c         hic(*)  - HIstory Correspondence table

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_geom.h'
      include  'c_keyh.h'
      include  'c_mate.h'
      include  'c_pair.h'
      include  'c_tole.h'

      integer   npair,hic((c_lp1+c_lp3),*), kv
      real*8    cs0(nr0,n0c1:nc01,*),cm0(nr0,n0c2:nc02,*)
      real*8    cp0(nr0,n0c3:nc03,*)

      save

      call cdebug0 ('    setcomp',-1)

c     Load control flags on common

      ifon   = cp0(3,1,npair)
      iffric = cp0(4,1,npair)
      ifsolm = cp0(2,2,npair)
      ifdeta = cp0(2,3,npair)
      ifaugm = cp0(2,5,npair)
      ifadhe = cp0(2,7,npair)

c     Load tolerances

      tlipen = cp0(3,6,npair)
      tlopen = cp0(4,6,npair)
      tlouts = cp0(5,6,npair)

c     Load offsets, set length of ch1,ch2,ch3, # of set, solution driver

      rnpair = cp0(1,-1,npair)
      ofh1   = cp0(2,-1,npair)
      ofh3   = cp0(3,-1,npair)
      lh1    = cp0(4,-1,npair)
      lh3    = cp0(5,-1,npair)
      nset   = cp0(6,-1,npair)
      ndrv   = abs(cp0(1, 0,npair))

c     Load # of element active till now

      nacte  = cp0(11,-1,npair)

c     Set contact search dimension

      cndm   = cp0(13,-1,npair)

c     Load surfaces information

      nsurf1 = abs(cp0(7,-1,npair))

      ofs1   = abs(cs0(2,-1,nsurf1))
      neps1  = abs(cs0(3,-1,nsurf1))
      dnope1 = abs(cs0(4,-1,nsurf1))
      ifsty1 = abs(cs0(1, 0,nsurf1))
      nope1  = abs(cs0(2, 0,nsurf1))

      nsurf2 = abs(cp0(8,-1,npair))

      ofs2   = abs(cs0(2,-1,nsurf2))
      neps2  = abs(cs0(3,-1,nsurf2))
      dnope2 = abs(cs0(4,-1,nsurf2))
      ifsty2 = abs(cs0(1, 0,nsurf2))
      nope2  = abs(cs0(2, 0,nsurf2))

c     Load material information

      nmat1  = cp0(9,-1,npair)
      if (nmat1.ne.0) then
        ifmty1 = cm0(1,0,nmat1)
        ofm1   = cm0(2,-1,nmat1)
      else
        ifmty1 = 0
        ofm1   = 0
      endif

c     Do the same only if mat # 2 is used

      nmat2  = cp0(10,-1,npair)
      if (nmat2.ne.0) then
        ifmty2 = cm0(1,0,nmat2)
        ofm2   = cm0(2,-1,nmat2)
      else
        ifmty2 = 0
        ofm2   = 0
      endif

c     Load correspondence table for history variables

      do kv = 1,c_lp1
        p1(kv) = hic(kv,npair)
      end do
      do kv = 1,c_lp3
         p3(kv) = hic(c_lp1+kv,npair)
      end do

      end
