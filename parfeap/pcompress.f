c$Id: pcompress.f,v 1.1 2006/11/21 16:44:39 rlt Exp $
      subroutine pcompress( ld, ls, p, s, nst )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate compressed stiffness/rhs
c
c      Inputs:
c         lr(*)    - Row numbers active
c         lc(*)    - Column numbers active
c         p(*)     - Full vector
c         s(nst,*) - Full matrix
c         nst      - Size of matrix dimensioned

c      Outputs:
c         p(*)     - Compressed vector
c         s(nst,*) - Compressed matrix
c         nar      - Number of active rows
c         nac      - Number of active columnts
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'pfeapb.h'

      integer    nst, m,n
      integer    ld(nst,*), ls(nst)
      real*8     p(nst,2), s(nst,nst,*)

c     Build row compression list

      nar = 0
      do n = 1,nst
        if(ld(n,1).gt.0) then
           nar   = nar + 1
           ls(nar)   = n
        endif
      end do ! n

c     Compress array

      if(nar.ne.nst) then
        do n = 1,nar
          ld(n,1) = ld(ls(n),1)
          ld(n,3) = ld(ls(n),3)
          ld(n,4) = ld(ls(n),4)
          ld(n,5) = ld(ls(n),5)
          ld(n,6) = ld(ls(n),6)
          p(n,1)    = p(ls(n),1)
        end do ! n
        do n = nar+1,nst
          ld(n,1) = 0
          ld(n,3) = 1
          ld(n,4) = 0
          ld(n,5) = 0
        end do
        do m = 1,nst
          do n = 1,nar
            s(n,m,1) = s(ls(n),m,1)
          end do ! n
        end do ! m
 
      endif

c     Build column compression list

      nac = 0
      do n = 1,nst
        if(ld(n,2).gt.0) then
          nac     = nac + 1
          ls(nac) = n
        endif
      end do ! n

c     Compress array

      if(nac.ne.nst) then
        do m = 1,nac
          ld(m,2) = ld(ls(m),2)
          do n = 1,nar
            s(n,m,1) = s(n,ls(m),1)
          end do ! n
        end do ! m
        do m = nac+1,nst
          ld(m,2) = 0
        end do ! m

      endif

c     Set compressed assembly

      do n = 1,nar
        ld(n,3) = 1
        ld(n,4) = nac
      end do ! n

      end ! pcompress
