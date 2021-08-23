c$Id: membr3d.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine membr3d(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Quadrilateral membrane element for feap

c     Input parameters set as follows:

c       Small and finite deformation
c         ndm = 3 (x,y,z cartesian coordinates at nodes)
c         ndf = 3 (u-x,u-y,u-z at nodes)
c         nen = 4 nodes (counterclockwise around element)

c       Note: 1-direction bisects diagonals between 2-3 element and
c             2-direction bisects diagonals between 3-4 element nodes.
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      implicit  none

      include  'eldata.h'
      include  'hdata.h'
      include  'iofile.h'
      include  'mdata.h'
      include  'strnum.h'

      include  'comblk.h'

      integer   ndm,ndf,nst,isw, tdof, i,j, ix(*)

      real*8    d(*),xl(ndm,*),ul(ndf,*),s(nst,*),p(nst)

      save

c     Input material properties

      if(isw.eq.1) then

        if(ior.lt.0) write(*,2000)
        write(iow,2000)
        call inmate(d,tdof, 12 ,5)

c       Deactivate dof in element for dof > 3

        do i = 4,ndf
          ix(i) = 0
        end do ! i

c       Set plot sequence

        pstyp = 2

c       Initialize the finite deformation shell

        if(d(18).lt.0.0d0) then
          call mem3df(d,ul,xl,s,p,ndf,ndm,nst,isw)
        endif

c       Set maximum plot variable number

        istv = 15

c     Initialize element strains for activation

      elseif(isw.eq.17) then

        j = nh3
        do i = 1,nel
          hr(j  ) = ul(1,i)
          hr(j+1) = ul(2,i)
          hr(j+2) = ul(3,i)
          j     = j + 3
        end do ! i

c     Initialize element strains for deactivation

      elseif(isw.eq.18) then

        do i = 1,3*nel
          hr(nh3+i-1) = 0.0d0
        end do ! i

c     External node check

      elseif(isw.eq.26) then

        call pcorner2d()

C     Remaining options

      else

c       Small deformation

        if(d(18).gt.0.0d0) then
          call mem3ds(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)

c       Finite deformation

        else
          call mem3df(d,ul,xl,s,p,ndf,ndm,nst,isw)
        endif

      endif

c     Format

2000  format(5x,'T h r e e   D i m e n s i o n a l   M e m b r a n e',
     &          '   E l e m e n t')

      end
