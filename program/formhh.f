c$Id: formhh.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine formhh(hh,g,neqg,neq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Forms finite element G and H arrays

c      Inputs:
c         neqg   - Number of constraint equations
c         neq    - Number of active equations

c      Outputs:
c         hh(neqg,*) - Export matrix  (TEMP5: hr(np(115)))
c         g(neq,*,2) - G-vector       (TEMP6: hr(np(116)))
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'eqsym.h'
      include  'part0.h'
      include  'pointer.h'
      include  'pscal.h'
      include  'comblk.h'

      integer   neqg,neq,mm,nn, i
      real*8    engy
      real*8    hh(neqg,neqg),g(neq,neqg,2)

      save

c     Loop through equations to form Gu' = A_inv * Gu

      do mm = 1,neqg

        call dasol(hr(np(4+npart)),hr(np(npart)+neq),
     &             hr(np(npart)),g(1,mm,1),mr(np(20+npart)),
     &             neqs,neq,engy,scale(npart))

        do nn = 1,neqg
          do i = 1,neq
            hh(nn,mm) = hh(nn,mm) - g(i,nn,2)*g(i,mm,1)
          end do ! i
        end do ! nn
      end do ! mm

      end
