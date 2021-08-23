c$Id: getincre.f,v 1.1 2006/11/20 20:32:46 rlt Exp $
      subroutine getincrem (ch2,ndats,kset,kc,gnv,dgn,dfn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996           1.00
c               Robert L. Taylor         October 12, 1996           1.01

c      Acronym: AUGMENTation

c      Purpose: Augment contact forces

c      Inputs :
c         ch2(*)  - Contact history variables (current)

c      Outputs:
c         ch2(*)  - Contact history variables (current)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_contac.h'
      include  'c_keyh.h'

      integer   ndats,kset,kc, istgn,istgn1,istgn2,ks,kr
      real*8    ch2(*), dgn(ndats,*),dfn(ndats,*),gnv(*)
      real*8    fn,gn,fn1,gn1,fn2,gn2, fndp,fn1dp, gndp,gn1dp

      save

      call cdebug0 ('      getincrem',-1)

c     data shift

      do kr = ndats,2,-1
        dgn(kr,kset) = dgn(kr-1,kset)
        dfn(kr,kset) = dfn(kr-1,kset)
      end do

c     check for opening

      istgn  = ch2(p1(4))

      if (istgn.ge.0) then
        fn  = ch2(p1(51))
        fn1 = ch2(p1(51)+1)
        fn2 = ch2(p1(51)+2)

        gn  = ch2(p1(9))
        gn1  = ch2(p1(9)+1)
        gn2  = ch2(p1(9)+2)

        rewind (77)
        write (77,*) gn
        rewind (77)
        read (77,*) gndp

        rewind (77)
        write (77,*) gn1
        rewind (77)
        read (77,*) gn1dp

        rewind (77)
        write (77,*) fn
        rewind (77)
        read (77,*) fndp

        rewind (77)
        write (77,*) fn1
        rewind (77)
        read (77,*) fn1dp

c       Store

        if (fn1.ne.0.d0) then
          gnv(kset)   = gn
          dgn(1,kset) = gn - gn1
          dfn(1,kset) = fn - fn1
        else
          gnv(kset)   = 0.d0
          dgn(1,kset) = 0.d0
          dfn(1,kset) = 0.d0
        endif

c       Update values on element history set

        ch2(p1(51)+1) = fn
        ch2(p1(51)+2) = fn1
        ch2(p1(9)+1)  = gn
        ch2(p1(9)+2)  = gn1

c     Deal with opening

      else
        gnv(kset)   = 0.d0
        dgn(1,kset) = 0.d0
        dfn(1,kset) = 0.d0

        ch2(p1(51)+1) = 0.d0
        ch2(p1(51)+2) = 0.d0
        ch2(p1(9)+1)  = 0.d0
        ch2(p1(9)+2)  = 0.d0
      endif

c     Update normal status flag

      istgn1 = ch2(p1(4)+1)
      istgn2 = ch2(p1(4)+2)

      ch2(p1(4)+1) = istgn
      ch2(p1(4)+2) = istgn1

      if (ifdb) then
        if (kset.eq.nset) then
          write (69,*) 'DGN, column ',kc
          do ks = 1,nset
            write (69,*) dgn(kc,ks)
          end do
          write (69,*) 'DFN'
          do ks = 1,nset
            write (69,*) dfn(kc,ks)
          end do
          write (69,*) 'GN'
          do ks = 1,nset
            write (69,*) gnv(ks)
          end do
        endif
      endif

      end
