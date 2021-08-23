c$Id: ptlocal.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine ptlocal(ul,p,s, cplxfl, ndf,nen,nel,nst, nrot, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Rotate element solution parameters and arrays for inclined
c               boundary nodes

c      Inputs:
c        ul(*)  - Solution parameters
c        p(*)   - Element vector
c        s(*)   - Element array
c        cplxfl - Complex flag
c        nrot(2)- Number of rotations: (1) Angle; (2) Euler angles
c        ndf    - Dimension 1 of ul
c        nen    - Dimension 2 of ul
c        nel    - Number of nodes on element
c        nst    - Size of stiffness
c        isw    - Switch parameter

c      Outputs:
c        ul(*)  - Rotated solution parameters
c        p(*)   - Rotated element vector
c        s(*)   - Rotated element array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'mdata.h'
      include   'pointer.h'
      include   'comblk.h'

      logical    cplxfl
      integer    nrot(2), ndf,nen,nel,nst, isw
      real*8     ul(ndf,nen,*),p(*),s(nst,*)

      save

c     Planar angle rotations

      if(nrot(1).gt.0) then
        call ptrans(dal,hr(np(46)),ul,p,s,nel,ndf,nst,isw)
        if(ral(1).ne.0) then
          call ptrans(ral,hr(np(46)),ul,p,s,nel,ndf,nst,isw)
        endif
        if(cplxfl) then
          call ptrans(dal,hr(np(46)),ul(1,1,8),p(nst+1),
     &                s(nst,nst+1),nel,ndf,nst,isw)
          if(ral(1).ne.0) then
            call ptrans(ral,hr(np(46)),ul(1,1,8),p(nst+1),
     &                s(nst,nst+1),nel,ndf,nst,isw)
          endif
        endif
      endif

c     Euler angle transformations

      if(nrot(2).gt.0) then
        call petrans(dal,hr(np(243)),ul,p,s,nel,ndf,nst,isw)
        if(ral(1).ne.0) then
          call petrans(dal,hr(np(243)),ul,p,s,nel,ndf,nst,isw)
        endif
        if(cplxfl) then
          call petrans(dal,hr(np(243)),ul(1,1,8),p(nst+1),
     &                 s(nst,nst+1),nel,ndf,nst,isw)
          if(ral(1).ne.0) then
            call petrans(dal,hr(np(243)),ul(1,1,8),p(nst+1),
     &                 s(nst,nst+1),nel,ndf,nst,isw)
          endif
        endif
      endif

      end
