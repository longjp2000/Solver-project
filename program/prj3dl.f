c$Id: prj3dl.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine prj3dl(pl,ix,x,x0,id,ie,f,ibn,ang,gap0,nen,nen1,ndm,
     &                  ndf,numnp,numel,prt,prth,fnorm,polfl,ddof,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use 'pi' from 'pconstant.h'                        14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute loads on 3-D surfaces

c      Inputs:
c         pl(ndf)        - Patterns
c         ix(nen1,*)     - Element nodal connection lists
c         ie(*)          - Element information array
c         x(ndm,*)       - Nodal coordinates
c         x0(ndm)        - Origin of polar coordinates
c         gap0           - User specified gap
c         nen            - Maximum number of nodes on any element
c         nen1           - Dimension of ix array
c         ndm            - Spatial dimension of mesh
c         ndf            - Degree of freedoms/node
c         numnp          - Number of nodes
c         numel          - Number of elements/faces
c         prt            - Output values when true
c         prth           - Output header when true
c         fnorm          - Force/Displ type: 1 = Normal,
c                                            2 = Tangential,
c                                            3 = Displacement
c         polfl          - Polar/cylindrical coordinate flag
c         ddof           - Direction for displacement/tangent
c         isw            - Switch:           1 = Force/Displ.;
c                                            2 = Boundary codes
c                                            3 = Angle values

c      Temporary:
c         ibn(numnp)     - Marks nodes on patch

c      Outputs:
c         f(ndf,numnp,*) - Nodal loads          : isw = 1
c         id(ndf,*)      - Nodal boundary codes : isw = 2
c         ang(*)         - Angle boundary cond  : isw = 3
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdat1.h'
      include  'iofile.h'
      include  'pconstant.h'
      include  'pdata5.h'
      include  'pdata6.h'
      include  'trdata.h'

      logical   polfl, errc, pinput, prt, prth
      integer   i,j,k,l,n, i1,i2, iel,iiel,ifac,lint,projpt,ma,nel,pstyp
      integer   nef,nev,nen,nen1,ndm,ndf,numnp,numel,fnorm,ddof,isw
      integer   ibn(numnp), ix(nen1,numel), id(ndf,numnp), ie(nie,*)
      integer   ic(9), iv(4,7), iq(5,7), it(3,4), il(9,7), itq(6,4)
      real*8    pl(*),x(ndm,numnp),x0(ndm),f(ndf,numnp,*),ang(*)
      real*8    gap0,gap,tol,tolgap,load,xsj(3),dir,den,th
      real*8    td(5),xmin(3),xmax(3),xp(3,9),v1(3),v2(3),xi(3),xx(3)
      real*8    nn(3),normal(3),pr(9),pp(9),shp(9),shps(9),xl(3,9)
      real*8    sg(3,9),tg(4,7),fl(3,9),xc(3),vv(2), angl,sn,cn

      save

c     8-node brick faces

      data iv  /1,4,3,2, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 5,6,7,8,
     &          1,2,3,4/

c     20 and 27-node brick sides and center

      data iq  /12,11,10, 9,21,  9,18,13,17,25, 10,19,14,18,24,
     &          11,20,15,19,26, 12,17,16,20,23, 13,14,15,16,22,
     &           9,10,11,12,21/

c     4-node tetrahedron faces

      data it  /1,2,4, 2,3,4, 3,1,4, 1,3,2/

c     10-node tetrahedron faces

      data itq / 1, 2, 4, 5, 9, 8, 2, 3, 4, 6,10, 9,
     &           3, 1, 4, 7, 8,10, 1, 2, 3, 5, 6, 7 /
c     Solution tolerances

      data tol /1.0d-3/

c     Input points

      do n = 1,9
        ic(n) = 0
        pp(n) = 0.d0
      end do ! n

c     Force-pressure-displacements inputs

      if(isw.eq.1) then
        if(fnorm.eq.1) then          ! Normal pressure/force
          write(iow,2000) (n,n=1,3)
        elseif(fnorm.eq.3) then      ! Displacement component
          write(iow,2010) (n,n=1,3)
        endif
        i = 1
        do while(i.ne.0)
          errc = pinput(td,5)
          if(errc .or. td(1).eq.0.0d0) then
            backspace(ior)
            i       = 0
          else
            i       = td(1)
            xp(1,i) = td(2)
            xp(2,i) = td(3)
            xp(3,i) = td(4)
            pp(i)   = td(5)
            ic(i)   = i
            write(iow,2001) i,(xp(n,i),n=1,3),pp(i)
          endif
        end do ! i.ne.0

c     Boundary Code inputs

      elseif(isw.eq.2) then

        write(iow,2004) (n,n=1,3),(n,n=1,ndf)
        write(iow,2005) (int(pl(n)),n=1,ndf)
        i = 1
        do while(i.ne.0)
          errc = pinput(td,5)
          if(errc .or. td(1).eq.0.0d0) then
            backspace(ior)
            i       = 0
          else
            i       = td(1)
            xp(1,i) = td(2)
            xp(2,i) = td(3)
            xp(3,i) = td(4)

            ic(i)   = i
            write(iow,2001) i,(xp(n,i),n=1,3)
          endif
        end do ! i.ne.0

c     Boundary angle inputs

      elseif(isw.eq.3) then

        write(iow,2011) (n,n=1,3)
        i = 1
        do while(i.ne.0)
          errc = pinput(td,5)
          if(errc .or. td(1).eq.0.0d0) then
            backspace(ior)
            i       = 0
          else
            i       = td(1)
            xp(1,i) = td(2)
            xp(2,i) = td(3)
            xp(3,i) = td(4)
            pp(i)   = td(5)
            ic(i)   = i
            write(iow,2001) i,(xp(n,i),n=1,3),pp(i)
          endif
        end do ! i.ne.0
      end if

c     Set missing points

      if(min(ic(1),ic(2),ic(3),ic(4)).eq.0) then
        write(ilg,3000) ' -->Missing corner coordinate:'
        write(ilg,3000) '    Check surface node numbers'
        write(iow,3000) ' -->Missing corner coordinate:'
        write(iow,3000) '    Check surface node numbers'
        call plstop()
      endif

c     Fill missing mid-side nodes

      do i = 1,4
        if(ic(i+4).eq.0) then
          j         = mod(i,4) + 1
          xp(1,i+4) = 0.5d0*(xp(1,i) + xp(1,j))
          xp(2,i+4) = 0.5d0*(xp(2,i) + xp(2,j))
          xp(3,i+4) = 0.5d0*(xp(3,i) + xp(3,j))
          pp(i+4)   = 0.5d0*(  pp(i) +   pp(j))
          j = i+4
          if(isw.eq.1 .or. isw.eq.3) then
            write(iow,2001) j,(xp(n,j),n=1,3),pp(j)
          elseif(isw.eq.2) then
            write(iow,2001) j,(xp(n,j),n=1,3)
          endif
        endif
      end do ! i

c     Fill central node

      if(ic(9).eq.0) then
        xp(1,9) = 0.50d0*(xp(1,5) + xp(1,6) + xp(1,7) + xp(1,8))
     &          - 0.25d0*(xp(1,1) + xp(1,2) + xp(1,3) + xp(1,4))
        xp(2,9) = 0.50d0*(xp(2,5) + xp(2,6) + xp(2,7) + xp(2,8))
     &          - 0.25d0*(xp(2,1) + xp(2,2) + xp(2,3) + xp(2,4))
        xp(3,9) = 0.50d0*(xp(3,5) + xp(3,6) + xp(3,7) + xp(3,8))
     &          - 0.25d0*(xp(3,1) + xp(3,2) + xp(3,3) + xp(3,4))
        pp(9)   = 0.50d0*(pp(5)   + pp(6)   + pp(7)   + pp(8)  )
     &          - 0.25d0*(pp(1)   + pp(2)   + pp(3)   + pp(4)  )
        i = 9
        if(isw.eq.1 .or. isw.eq.3) then
          write(iow,2001) i,(xp(n,i),n=1,3),pp(i)
        elseif(isw.eq.2) then
          write(iow,2001) i,(xp(n,i),n=1,3)
        endif
      endif

c     Transform surface coordinates to Global frame

      do n = 1,9
        xc(1)   = xp(1,n)
        xc(2)   = xp(2,n)
        xc(3)   = xp(3,n)
        xp(1,n) = xr(1) + tr(1,1)*xc(1) + tr(1,2)*xc(2) + tr(1,3)*xc(3)
        xp(2,n) = xr(2) + tr(2,1)*xc(1) + tr(2,2)*xc(2) + tr(2,3)*xc(3)
        xp(3,n) = xr(3) + tr(3,1)*xc(1) + tr(3,2)*xc(2) + tr(3,3)*xc(3)
      end do ! n

c     Initialize surface node indicators to zero

      do n = 1,numnp
        ibn(n) = 0
      end do ! n

c     Set min/max for patch coordinates

      do i = 1,3
        xmin(i) = xp(i,1)
        xmax(i) = xp(i,1)
        do n = 2,9
          xmin(i) = min(xmin(i),xp(i,n))
          xmax(i) = max(xmax(i),xp(i,n))
        end do ! n
      end do ! i

c     Polar angle conversion

      th = 180.0d0/pi

c     Set gap value

      gap = 0.0d0
      do i = 1,3
        gap = gap + xmax(i) - xmin(i)
      end do ! i
      gap = 1.d-4*gap

      if(gap0 .gt. 0.0d0) then
        gap    = gap0
      else
        tolgap = tol/dble(numnp)**(1.d0/3.d0)
        do i = 1,3
          gap = max(gap,tolgap*(xmax(i)-xmin(i)))
        end do ! i
      endif

      do i = 1,3
        xmin(i) = xmin(i) - 1.5d0*gap
        xmax(i) = xmax(i) + 1.5d0*gap
      end do ! i

c     Determine potential nodes for loads, etc.

      if(polfl) then
        angl  = 0.5d0*(xmin(2) + xmax(2))
        vv(1) = cos(angl/th)
        vv(2) = sin(angl/th)
      endif
      do n = 1,numnp
        if(polfl) then
          xc(1) =  sqrt((x(1,n) - x0(1))**2 + (x(2,n) - x0(2))**2)
          cn    =  vv(1)*(x(1,n) - x0(1)) + vv(2)*(x(2,n) - x0(2))
          sn    = -vv(2)*(x(1,n) - x0(1)) + vv(1)*(x(2,n) - x0(2))
          xc(2) =  atan2(sn,cn)*th + angl
          xc(3) =  x(3,n) - x0(3)
        else
          xc(1) =  x(1,n)
          xc(2) =  x(2,n)
          xc(3) =  x(3,n)
        endif

c       Eliminate nodes outside patch block

        if((xc(1).lt.xmin(1) .or. xc(1).gt.xmax(1)) .or.
     &     (xc(2).lt.xmin(2) .or. xc(2).gt.xmax(2)) .or.
     &     (xc(3).lt.xmin(3) .or. xc(3).gt.xmax(3)) ) then
          ibn(n) = -1
        endif

c       Search remaining nodes for those near patch

        if(ibn(n).ge.0) then
          ibn(n) = projpt(xc,xp,xi,gap,normal,shp)
        endif
      end do ! n

c     Force computations for normal pressures

      if(isw.eq.1) then

c       Determine elements with faces on patch

        if(prt) then
          call prtitl(prth)
          if(fnorm.eq.1) then
            write(iow,2002) (j,j=1,3)
          elseif(fnorm.eq.3) then
            write(iow,2008)
          endif
        endif

        do n = 1,numel

c         Determine element type

          ma    = ix(nen1,n)
          pstyp = ie(1,ma)

          if(pstyp.gt.0) then
            iel   = ie(nie-1,ma)
            do i = nen,1,-1
              if(ix(i,n).gt.0) then
                nel = i
                exit
              endif
            end do ! i

            call plftyp(pstyp,nel,iel)
            if(iel.ge.0) then
              iiel = inord(iel)
            else
              iiel = exord(-iel)
            endif

c           No face if iiel < 0

            if(iiel.lt.0) then

c           Linear tetrahedral element faces

            elseif(iiel .eq. 9 ) then
  
              nef = 3
              nev = 3
              i1  = 4
              i2  = 1
              do i = 1,4
                do j = 1,3
                  il(j,i) = it(j,i)
                end do ! j
              end do ! i

c           Quadratic tetrahedral element faces

            elseif(iiel .eq. 15 ) then

              nef = 6
              nev = 3
              i1  = 4
              i2  = 1
              do i = 1,4
                do j = 1,6
                  il(j,i) = itq(j,i)
                end do ! j
              end do ! i

c           Brick and shell element faces

            else

              i1  = 6
              i2  = 1
              nev = 4
              if(nel.eq.27) then
                nef = 9
              elseif(nel.eq.20) then
                nef = 8
              elseif(nel.eq.8) then
                nef = 4
              elseif(nel.eq.4) then
                nef = 4
                i1  = 7
                i2  = 6
              else
                nef = 4
                i1  = 0
                i2  = 1
              endif
              do i = 1,7
                do j = 1,4
                  il(j,i) = iv(j,i)
                end do ! j
                do j = 5,nef
                  il(j,i) = iq(j-4,i)
                end do ! j
              end do ! i
  
            endif

            do i = 1,i1,i2
              ifac = 0
              do j = 1,nef
                ic(j) = ix(il(j,i),n)
                if(ic(j).gt.0) then
                  ifac = ifac + ibn(ic(j))
                endif
              end do ! j

c             Compute pressures for matching faces

              if(ifac.eq.nef) then
                do j = 1,3
                  v1(j) =  x(j,ic(2)) - x(j,ic(nev))
                  v2(j) =  x(j,ic(3)) - x(j,ic(1))
                  xx(j) =  0.0d0
                  do k = 1,nef
                    xx(j)   = xx(j) + x(j,ic(k))
                    xl(j,k) = x(j,ic(k))
                  end do ! k
                  xx(j) = xx(j)/dble(nef)
                end do ! j
                nn(1) = v1(2)*v2(3) - v1(3)*v2(2)
                nn(2) = v1(3)*v2(1) - v1(1)*v2(3)
                nn(3) = v1(1)*v2(2) - v1(2)*v2(1)
                if(polfl) then
                  xc(1) =  sqrt((xx(1) - x0(1))**2 + (xx(2) - x0(2))**2)
                  cn    =  vv(1)*(xx(1) - x0(1)) + vv(2)*(xx(2) - x0(2))
                  sn    = -vv(2)*(xx(1) - x0(1)) + vv(1)*(xx(2) - x0(2))
                  xx(2) =  atan2(sn,cn)*th + angl
                  xx(1) =  xc(1)
                  call pdegree(xx(2), sn,cn)
                  xx(3) =  xx(3) - x0(3)
                  xc(1) =  nn(1)*cn + nn(2)*sn
                  nn(2) = -nn(1)*sn + nn(2)*cn
                  nn(1) =  xc(1)
                endif
                ifac  = projpt(xx,xp,xi,gap,normal,shp)
                dir   = nn(1)*normal(1)+nn(2)*normal(2)+nn(3)*normal(3)
                den   = sqrt(nn(1)*nn(1)+nn(2)*nn(2)+nn(3)*nn(3))
     &                * sqrt(normal(1)*normal(1)+normal(2)*normal(2)
     &                      +normal(3)*normal(3))
                if(dir.gt.0.8d0*den) then

c                 Compute load at element nodes

                  do j = 1,nef
                    if(polfl) then
                      xc(1) =  sqrt((xl(1,j) - x0(1))**2
     &                            + (xl(2,j) - x0(2))**2)
                      cn    =  vv(1)*(x(1,j) - x0(1))
     &                       + vv(2)*(x(2,j) - x0(2))
                      sn    = -vv(2)*(x(1,j) - x0(1))
     &                       + vv(1)*(x(2,j) - x0(2))
                      xc(2) =  atan2(sn,cn)*th + angl
                      xc(3) =  xl(3,j) - x0(3)
                    else
                      xc(1) = xl(1,j)
                      xc(2) = xl(2,j)
                      xc(3) = xl(3,j)
                    endif
                    ifac  = projpt(xc,xp,xi,gap,normal,shp)
                    pr(j) = 0.0d0
                    do k = 1,9
                      pr(j) = pr(j) + shp(k)*pp(k)
                    end do ! k
                    fl(1,j) = 0.0d0
                    fl(2,j) = 0.0d0
                    fl(3,j) = 0.0d0
                  end do ! j

c                 Normal Force:

                  if(fnorm.eq.1) then

c                   Compute loads on nodes

                    if(nef.eq.4) then
                      call int2d(2,lint,sg)
                    elseif(nef.ge.8) then
                      call int2d(3,lint,sg)
                    else
                      if(nef.eq.3) then
                        l = 3
                      elseif(nef.eq.6) then
                        l = 7
                      endif
                      call tint2d(l,lint,tg)
                      do l = 1,lint
                        sg(1,l) = tg(1,l)
                        sg(2,l) = tg(2,l)
                        sg(3,l) = tg(4,l)*0.5d0
                      end do ! l
                    endif
                    do l = 1,lint
  
                      call shp3p(sg(1,l),xl,shps,xsj,nef)
                      load  = 0.0d0
                      do k = 1,nef
                        load  = load  + shps(k)*pr(k)
                      end do ! k
  
                      load = load*sg(3,l)
                      do k = 1,nef
                        fl(1,k) = fl(1,k) + load*xsj(1)*shps(k)
                        fl(2,k) = fl(2,k) + load*xsj(2)*shps(k)
                        fl(3,k) = fl(3,k) + load*xsj(3)*shps(k)
                      end do ! k
                    end do ! l

c                   Check for sloping boundaries

                    do k = 1,nef
                      if(ang(ic(k)).ne.0.0d0) then
                        call pdegree(ang(ic(k)), sn,cn)
                        td(3)   =  cn*fl(1,k) + sn*fl(2,k)
                        fl(2,k) = -sn*fl(1,k) + cn*fl(2,k)
                        fl(1,k) =  td(3)
                      endif
                    end do ! k

                    if(prt) then
                      do k = 1,nef
                        write(iow,2003) k,ic(k),(fl(j,k),j=1,3),pr(k)
                      end do ! k
                    endif

                    do k = 1,nef
                      f(1,ic(k),1) = f(1,ic(k),1) + fl(1,k)
                      f(2,ic(k),1) = f(2,ic(k),1) + fl(2,k)
                      f(3,ic(k),1) = f(3,ic(k),1) + fl(3,k)
                    end do ! k

c                 Displacement: Component ddof

                  elseif(fnorm.eq.3) then

                    do k = 1,nef
                      f(ddof,ic(k),2) = pr(k)
                      if(prt) then
                        write(iow,2009) k,ic(k),ddof,f(ddof,ic(k),2)
                        if(ior.lt.0) then
                          write(*,2009) k,ic(k),ddof,f(ddof,ic(k),2)
                        endif
                      endif
                    end do ! k

                  end if ! fnorm

                endif ! dir
              endif ! ifac

            end do ! i
          endif ! pstyp > 0
        end do ! n

c     Set the boundary codes

      elseif(isw.eq.2) then

        if(prt) then
          call prtitl(prth)
          write(iow,2006) (j,j=1,ndf)
        endif
        do n = 1,numnp
          if(ibn(n).eq.1) then
            do j = 1,ndf
              id(j,n) = abs(id(j,n)) + int(pl(j))
            end do ! j
            if(prt) then
              write(iow,2007) n,(id(j,n),j=1,ndf)
            endif
          endif
        end do ! n

c     Set the boundary angles

      elseif(isw.eq.3) then
        if(prt) then
          write(iow,2012)
        endif
        do n = 1,numnp
          if(ibn(n).eq.1) then

c           Interpolate angle on patch to node position

            ifac  = projpt(xc,xp,xi,gap,normal,shp)
            pr(1) = 0.0d0
            do j = 1,9
              pr(1) = pr(1) + shp(j)*pp(j)
            end do ! j
            if(prt) then
              write(iow,2013) n,pr(1)
            endif

c           Store global angle

            ang(n) = pr(1)

          endif
        end do ! n

      endif

c     Formats

2000  format('   C o o r d i n a t e    S u r f a c e   L o a d s'/
     &       /6x,'Node',3(i5,' Coord'),3x,'Pressure')
2001  format(i10,1p,6e11.3)
2002  format(/7x,'N o d a l    F o r c e s'//3x,'Local',5x,'Global'/
     &       4x,'Node',6x,'Node',3(i5,' Force'),3x,'Pressure')
2003  format(i8,i10,1p,5e11.3)

2004  format(/6x,'Node',3(i5,' Coord'),9(i2,' BC':))
2005  format(43x,9i5)

2006  format(4x,'Node',9(i3,' BC':))
2007  format(i10,9i6)

2008  format(/7x,'N o d a l    D i s p l a c e m e n t s'
     &      //3x,'Local',3x,'Global'/
     &        4x,'Node',4x,'Node     DOF Displacement')

2009  format(3i8,1p,e13.4)

2010  format('   C o o r d i n a t e    S u r f a c e',
     &       '   D i s p l a c e m e n t s'/
     &       /6x,'Node',3(i5,' Coord'),5x,'Displ.')

2011  format('   C o o r d i n a t e    S u r f a c e   A n g l e s'/
     &       /6x,'Node',3(i5,' Coord'),6x,'Angle')
2012  format(/7x,'N o d a l    A n g l e s'//
     &       10x,'Global Node    Angle'/)
2013  format(i20,1p,5e11.3)

3000  format(/' *ERROR* PRJ3DL: ',a)

      end
