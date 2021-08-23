c$Id: scalev.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine scalev(v,pdf,ndm,ndf,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Scale vector to have maximum element of +1.0

c      Inputs:
c         v(ndf,*) - Vector of values
c         pdf(*)   - DOF to scale on
c         ndm      - Space dimension of mesh
c         ndf      - DOF's/node (maximum)
c         numnp    - Number of nodes

c      Outputs:
c         v(ndf,*) - Unit vector of values
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,n,ndm,ndf,numnp, pdf(*)
      real*8    v(ndf,*),vmax

      save

c     Locate maximum

      vmax = 0.0d0
      do i = 1,ndm
        if(pdf(i).ge.1.and.pdf(i).le.ndf) then
          do n = 1,numnp
            vmax = max(vmax,abs(v(pdf(i),n)))
          end do ! n
        endif
      end do ! i

c     Perform scaling

      if(vmax.gt.0.0d0) then
        vmax = 1.d0/vmax
        do n = 1,numnp
          do i = 1,ndf
            v(i,n) = v(i,n)*vmax
          end do ! i
        end do ! n
      else
        write(*,*) 'Zero length vector in SCALEV'
      endif

      end
