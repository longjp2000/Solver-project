c$Id: updatl.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine updatl(id,ixt,du,u,ul,f,ndf,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Control program for Euler angle update 

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'corset.h'
      include   'counts.h'
      include   'iofile.h'
      include   'mdata.h'
      include   'part0.h'
      include   'ptest.h'
      include   'tdatb.h'

      include   'pointer.h'
      include   'comblk.h'

      integer    ndf,numnp, i,j,n
      real*8     cc
      integer    id(ndf,numnp,*),ixt(*)
      real*8     du(*),u(ndf,numnp,*),ul(ndf,3),f(ndf,numnp)
      real*8     ubl(20),ang

c     If test solution step set 'cc' to small value 'testva'

      if( testfl .and. niter.le.1 ) then
        cc = testva
      else
        cc = 1.d0 ! full step
      endif

      do n = 1,numnp
        do i = 1,ndf
          ubl(i) = u(i,n,1)
        end do ! i

c       Check for sloping angle boundaries

        if(anglefl) then
          if(hr(np(45)+n-1).ne.0.0d0) then
            ang = hr(np(45)+n-1)
            call upang(dal,ang,ubl,ndf,1)
            if(ndf.ge.6) then
              call upang(ral,ang,ubl,ndf,1)
            endif
          endif
        endif

c       Check for Euler angle boundaries

        if(eulerfl) then
          if(hr(np(242)+3*n-3).ne.0.0d0 .or.
     &       hr(np(242)+3*n-2).ne.0.0d0 .or.
     &       hr(np(242)+3*n-1).ne.0.0d0) then
            call upeul(dal,hr(np(242)+3*n-3),ubl,ndf,1)
            if(ndf.ge.6) then
              call upeul(ral,hr(np(242)+3*n-3),ubl,ndf,1)
            endif
          endif
        endif

c       Extract local displacements

        do i = 1,ndf
          if(ndfp(i).eq.npart) then
            j = id(i,n,1)
            if (j.gt.0) then

c             For active dof compute values from solution
c             where 'du(j)' is increment of 'u' for active dof 'j'.

              ul(i,1) = du(j)*cc1
              ul(i,2) = du(j)*cc2
              ul(i,3) = du(j)

c           Set for restrained dof (not a rigid body)

            elseif(ixt(n).le.0 .and. id(i,n,2).gt.0) then

c             Set value from input for restrained dof

c             ul(i,1) = (f(i,n) - u(i,n,1))*cc
              ul(i,1) = (f(i,n) - ubl(i))*cc

c             Incremental boundary solution

              if(cc3.ne.0.0d0) then
                ul(i,3) = ul(i,1)*cc3
              else
                ul(i,3) = 0.0d0
                if(ul(i,1).ne.0.0d0) then
                  write(iow,*)' WARNING - infinite acceleration'
                endif
              endif
              ul(i,2) = ul(i,3)*cc2
            else
              ul(i,1) = 0.0d0
              ul(i,2) = 0.0d0
              ul(i,3) = 0.0d0
            endif
          endif
        end do ! i

c       Check for sloping angle boundaries

        if(anglefl) then
          if(hr(np(45)+n-1).ne.0.0d0) then
            call upang(dal,hr(np(45)+n-1),ul,ndf,2)
            if(ndf.ge.6) then
              call upang(ral,hr(np(45)+n-1),ul,ndf,2)
            endif
          endif
        endif

c       Check for Euler angle boundaries

        if(eulerfl) then
          if(hr(np(242)+3*n-3).ne.0.0d0 .or.
     &       hr(np(242)+3*n-2).ne.0.0d0 .or.
     &       hr(np(242)+3*n-1).ne.0.0d0) then
            call upeul(dal,hr(np(242)+3*n-3),ul,ndf,2)
            if(ndf.ge.6) then
              call upeul(ral,hr(np(242)+3*n-3),ul,ndf,2)
            endif
          endif
        endif

c       Perform update
        do i = 1,ndf
          if(ndfp(i).eq.npart) then
            u(i,n,1) = u(i,n,1) + ul(i,1)
            u(i,n,2) = u(i,n,2) + ul(i,2)
            u(i,n,3) =            ul(i,3)
          endif
        end do ! i
      end do ! n

      end
