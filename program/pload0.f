c$Id: pload0.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pload0(f,f0,u,nn,p)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Form nodal load vector for current time and partition

c      Inputs:
c         f(*,2)    - Current load/displacement vector
c         u(*)      - Current solution state
c         nn        - Number of components in vectors
c         p         - Current total proportional load level

c      Outputs:
c         f0(*,1)   - Fixed solution level for subsequent solutions
c         f0(*,2)   - Set to current solution state
c         f(*,*)    - Set to zero
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'part0.h'
      include  'pointer.h'
      include  'prld1.h'
      include  'sdata.h'
      include  'comblk.h'

      integer   i,n,nn, ipro
      real*8    p, f(nn,*),f0(nn,*),u(nn)

      save

c     Set nodal forces and displacements for NEWF command

      do i = 1,ndf
        if(ndfp(i).eq.npart) then
          do n = i,nn,ndf

c           Set force values

            ipro = mr(np(29)+n-1)
            if(ipro.eq.0) then
              f0(n,1) = f(n,1)*p + f0(n,1)
            else
              f0(n,1) = f(n,1)*prldv(ipro) + f0(n,1)
            endif

c           Set displacement values

            f0(n,2) = u(n)

c           Zero variable loading component

            f(n,1)  = 0.0d0
            f(n,2)  = 0.0d0
          end do ! n
        endif
      end do ! i

      end
