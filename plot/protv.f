c$Id: protv.f,v 1.1 2006/11/20 20:33:12 rlt Exp $
      subroutine protv(nty,u,ndm,ndf,numnp, du)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute displacements in global Cartesian frame.

c      Inputs:
c         nty(*)    - Nodal type
c         u(ndf,*)  - Solution vector at nodes
c         ndm       - Spatial dimension of mesh
c         ndf       - Number dof/node
c         numnp     - Number of nodes in mesh

c      Outputs:
c         du(ndf,*) - Cartesian displacements at nodes
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ndm,ndf,numnp, i,n, nty(*)
      real*8    u(ndf,*),du(ndf,*)

      save

      do n = 1,numnp
        if(nty(n).ge.0) then
          do i = 1,ndf
            du(i,n) = u(i,n)
          end do ! i
        else
          do i = 1,ndf
            du(i,n) = 0.0d0
          end do ! i
        endif
      end do ! n

      end
