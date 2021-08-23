c$Id: sprojb.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine sprojb(b,v,t,h,neq,nv,imas)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute subspace projection of 'b' to form 'h'

c      Inputs:
c         b(*)     - Symmetric coefficient matrix for eigenproblem
c         v(neq,*) - Set of iteration vectors
c         neq      - Number of equations in B
c         nv       - Size of projected matrix
c         imas     - Mass type: 1 = consistent; 2 = diagonal.

c      Scratch:
c         t(neq)   - Working vector

c      Outputs:
c         h(*)     - Projected matrix V_trans * B * V
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pointer.h'
      include  'comblk.h'

      integer   neq,nv,imas, i,j,k
      real*8    dot,b(*),v(neq,*),t(*),h(*)

      save

c     Compute 'z' and 'b' projection to form 'h'

      do j = 1,nv

c       Consistent mass

        if(imas.eq.1) then
          call pzero(t,neq)
          call caprod(b(1),b(neq+1),v(1,j),t,mr(np(90)),mr(np(91)),neq)

c       Lumped mass

        else
          do i = 1,neq
            t(i) = v(i,j)*b(i)
          end do ! i
        endif

c       Project 'z' and 'v' vectors to form 'h'

        k = j*(j+1)/2
        do i = j,nv
          h(k) = dot(t,v(1,i),neq)
          k = k + i
        end do ! i
        do i = 1,neq
          v(i,j) = t(i)
        end do ! i
      end do ! j

      end
