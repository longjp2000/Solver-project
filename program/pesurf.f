c$Id: pesurf.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pesurf(id,ix,ip,ep,xin,x,f,ang,ndf,ndm,nen,nen1,
     &                  numnp,numel,prt,prth,jsw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use 'one3' from 'pconstant.h' instead of 'third'   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Input surface loads or displacements from specified
c               coordinates and pressures.  Surface b.c. from
c               specified coordinates and pattern

c      Inputs:
c         ix(nen1,*)  - Element nodal connection data
c         x(ndm,*)    - Nodal coordinates of mesh
c         ndf         - Number dof/node
c         ndm         - Spatial dimension of mesh
c         nen         - Maximum number of nodes/element
c         nen1        - Dimension of ix array
c         numnp       - Number of nodes in mesh
c         numel       - Number of elements in mesh
c         prt         - Output results if true
c         prth        - Output title/header data if true
c         jsw         - Switch controling type of data to generate
c                       1 = Surface load/displ values
c                       2 = Boundary conditions
c                       3 = Angle (sloping)
c                       4 = Displacement conditions
c                       5 = Force conditions
c                       6 = Proportional load number conditions
c                      10 = Base load conditions
c                      11 = Line load conditions

c      Scratch:
c         ip(*)       - Nodal integer list storage
c         ep(numel,*) - Element list storage
c         xin(*)      - Nodal real list storage

c      Outputs:
c         id(ndf,*)   - Nodal boundary restraint conditions
c         f(ndf,*)    - Nodal force and boundary values
c         ang(*)      - Angles for sloping boundary nodes
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'comfil.h'
      include  'conval.h'
      include  'cornum.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'ioincl.h'
      include  'pconstant.h'
      include  'pointer.h'
      include  'trdata.h'
      include  'comblk.h'

      include  'p_int.h'

      character ptype*15,fext*4, wd(12)*4, type*4
      logical   prt,prth,errck,pinput,tinput,pcomp,polfl, wdflg
      logical   setmem,palloc, oprt,oprth
      integer   wdlist,fnorm,ddof, ndf,ndm,nen,nen1,numnp,numel,isw,jsw
      integer   i,i1,i2,i3,j,m,n,nn,nel,ntyp,ntot,numprt,l,lint
      integer   jmin
      integer   id(ndf,numnp,*),ix(nen1,*),ip(*),ep(numel,*),nend(2,2)
      real*8    cn,sn, y1,y2, tol0,tol,tolxi, ff,df,xi, d,gap0(2),alpha
      real*8    ximin,ximax,xia,xib,xic,fa,fb,fxa,fxb,fya,fyb, pa,pb,pc
      real*8    xc,yc, x1,x2, xy, xi1,xi2, zeta, zeta1,zeta2
      real*8    xin(*),x(ndm,*),f(ndf,numnp,*),ang(*)
      real*8    xe(2,4),xl(3,3),dxl(2), xp(3),x0(3),pl(30)
      real*8    shp(2,4),sw(2,4), fl(2,4), td(5), xld(2),xqd(2)

      save

      data      tol0 / 1.0d-5/, tolxi / 1.0d-8 /

      data      wdlist /12/
      data      wd / 'gap ', 'node', 'line', 'quad', 'pola', 'cart',
     &               'disp', 'tang', 'norm', 'surf', 'segm', 'dofs'/

c     Coordinate surface loading, boundary condition, and angle inputs

      isw = jsw
      if     (isw.eq.1) then
        ntot = nsurf - 1
      else if(isw.eq.2) then
        ntot = nbouf - 1
      else if(isw.eq.3) then
        ntot = nangf - 1
      else if(isw.eq.4) then
        ntot = ndisf - 1
      else if(isw.eq.5) then
        ntot = nforf - 1
      else if(isw.eq.6) then
        ntot = nprof - 1
      else if(isw.eq.7) then
        ntot = ndamf - 1
      else if(isw.eq.8) then
        ntot = nmasf - 1
      else if(isw.eq.9) then
        ntot = nstif - 1
      else if(isw.eq.10) then
        ntot = nbasf - 1
      else if(isw.eq.11) then
        ntot = nforl - 1
      else
        ntot = -1
      end if

      do nn = 0,ntot

c       Default to normal loading

        fnorm  = 1
        numprt = 0

c       Set file extender and number of input items

        if(isw.eq.1) then
          fext = 'sl0'
          ntyp = max(3,ndm+ndf)
        elseif(isw.eq.2) then
          fext = 'bn0'
          ntyp = ndm+ndf
        elseif(isw.eq.3) then
          fext = 'an0'
          ntyp = ndm+1
        elseif(isw.eq.4) then
          fext = 'ds0'
          ntyp = ndm+ndf
        elseif(isw.eq.5) then
          fext = 'fr0'
          ntyp = ndm+ndf
        elseif(isw.eq.6) then
          fext = 'yp0'
          ntyp = ndm+ndf
        elseif(isw.eq.7) then
          fext = 'vd0'
          ntyp = ndm+ndf
        elseif(isw.eq.8) then
          fext = 'ms0'
          ntyp = ndm+ndf
        elseif(isw.eq.9) then
          fext = 'ks0'
          ntyp = ndm+ndf
        elseif(isw.eq.10) then
          fext = 'hb0'
          ntyp = ndm+ndf
        elseif(isw.eq.11) then
          fext = 'pf0'
          ntyp = ndm
        endif

c       Check for file numbers

        if(nn.le.9) then
          write(fext(3:3),'(i1)') nn
        elseif(nn.le.99) then
          write(fext(2:3),'(i2)') nn
        endif
        call pinpfl('PESURF',fext, type, 1)
        oprt    = prt
        oprth   = prth
        gap0(1) = 0.0d0
        gap0(2) = 1.0d0
        polfl   = .false.

c       Input type and pressure on surface

1       ptype = ' '
        errck = tinput(ptype,1,pl,ntyp)
        if(errck) then
          write(ilg,3003)
          write(iow,3003)
          call plstop()
        endif

c       Check for legal data

        wdflg = .false.
        do i = 1,wdlist
          if(pcomp(ptype,wd(i),4)) wdflg = .true.
        end do ! i

c       Delete file and go to next list

        if(.not.wdflg) then
          call pinpfl('PESURF',fext, type, 2)
          go to 300
        elseif( pcomp(ptype,'norm',4)) then
          fnorm = 1
          go to 1
        elseif( pcomp(ptype,'tang',4)) then
          fnorm = 2
          go to 1
        elseif( pcomp(ptype,'disp',4) .or. pcomp(ptype,'dofs',4)) then
          ddof  = max(1,min(ndf,nint(pl(1))))
          fnorm = 3
          go to 1
        elseif( pcomp(ptype,'node',4)) then
          if(isw.eq.1 .or. isw.eq.5) then
            call pfboun(pl,x,f(1,1,1),mr(np(190)),ndm,ndf,numnp,numprt,
     &                  type,prt,prth,'-force')
          elseif(isw.eq.2) then
            call pcboun(pl,x,id(1,1,2),mr(np(190)),ndm,ndf,numnp,numprt,
     &                  gap0,type,prt,prth,'-b.c. ')
          elseif(isw.eq.3) then
            call paboun(pl,x,ang,mr(np(190)),ndm,numnp,numprt,
     &                  prt,prth)
          elseif(isw.eq.4) then
            call pcdisp(pl,x,f(1,1,2),mr(np(190)),ndm,ndf,numnp,numprt,
     &                  gap0,prt,prth)
          elseif(isw.eq.6) then
            call pcprop(pl,x,mr(np(29)),mr(np(190)),ndm,ndf,numnp,
     &                  numprt,prt,prth)
          elseif(isw.eq.7) then
            call pfboun(pl,x,hr(np(86)),mr(np(190)),ndm,ndf,numnp,
     &                  numprt,type,prt,prth,'-damp ')
          elseif(isw.eq.8) then
            call pfboun(pl,x,hr(np(87)),mr(np(190)),ndm,ndf,numnp,
     &                  numprt,type,prt,prth,'-mass ')
          elseif(isw.eq.9) then
            call pfboun(pl,x,hr(np(88)),mr(np(190)),ndm,ndf,numnp,
     &                  numprt,type,prt,prth,'-stif ')
          elseif(isw.eq.10) then
            call pcboun(pl,x,mr(np(125)),mr(np(190)),ndm,ndf,numnp,
     &                  numprt,gap0,type,prt,prth,'-base ')
          endif
          go to 1
        elseif( pcomp(ptype,'cyli',4)) then
           call pcangle(pl,x,ang,ndm,numnp,numprt,prt,prth)
          go to 1
        elseif( pcomp(ptype,'surf',4)) then
          if(jsw.eq.4) then
            isw   = 1
            ddof  = max(1,min(ndf,nint(pl(1))))
            fnorm = 3
          elseif(jsw.eq.5) then
            isw   = 1
            if(fnorm.eq.3) then ! Error
              write(ilg,3005)
              write(iow,3005)
              call plstop()
            endif
          endif
          setmem = palloc(113,'TEMP3',numnp,1)
          call prj3dl(pl,ix,x,x0,id(1,1,2),mr(np(32)),f,mr(np(113)),
     &                hr(np(45)),gap0,nen,nen1,ndm,ndf,numnp,numel,
     &                prt,prth,fnorm,polfl,ddof,isw)
          setmem = palloc(113,'TEMP3',0,1)
          go to 1
        elseif( pcomp(ptype,'gap', 3)) then
          gap0(1) = pl(1)
          gap0(2) = pl(2)
          if(gap0(2).eq.0.0d0) then
            gap0(2) = 1.0d0
          endif
          if(prt) then
            write(iow,*) '  GAP =',gap0
            if(ior.lt.0) then
              write(*,*) '  GAP =',gap0
            endif
          endif
          go to 1
        elseif( pcomp(ptype,'pola',4)) then
          polfl = .true.
          do i = 1,3
            x0(i) = pl(i)
          end do ! i
          if(prt) then
            write(iow,2004) x0
            if(ior.lt.0) then
              write(*,2004) x0
            endif
          endif
          go to 1
        elseif( pcomp(ptype,'cart',4)) then
          polfl = .false.
          go to 1
        elseif( pcomp(ptype,'segm',4)) then
          if(jsw.eq.4) then
            isw   = 1
            ddof  = max(1,min(ndf,nint(pl(1))))
            fnorm = 3
          elseif(jsw.eq.5) then
            isw   = 1
            if(fnorm.eq.3) then ! Error
              write(ilg,3005)
              write(iow,3005)
              call plstop()
            endif
          endif
          i  = 1
          i1 = 0
          do while (i.ne.0)
            errck = pinput(td,5)
            if(errck) then
              backspace ior
              i = 0
            else
              i  = max(0,nint(td(1)))
              i1 = max(i,i1)
              if(i.gt.3) then
                write(ilg,3001) i
                write(iow,3001) i
                call plstop()
              elseif(i.gt.0) then
                xl(1,i) = td(2)
                xl(2,i) = td(3)
                if(ndm.eq.3) then
                  xl(3,i) = td(4)
                else
                  xl(3,i) = 0.0d0
                endif
                if(isw.ne.2 .and. isw.ne.6) pl(i)   = td(ndm+2)
              endif
            endif
          end do ! while

c         For linear segments average to get mid value/point

          if(i1.eq.2) then
            xl(1,3) = 0.5d0*(xl(1,1) + xl(1,2))
            xl(2,3) = 0.5d0*(xl(2,1) + xl(2,2))
            xl(3,3) = 0.5d0*(xl(3,1) + xl(3,2))
            if(isw.ne.2 .and. isw.ne.6) pl(3)   = 0.5d0*(pl(1) + pl(2))
          endif

c       Quadratic segment

        elseif( pcomp(ptype,'quad',4)) then
          if(jsw.eq.4) then
            isw   = 1
            ddof  = max(1,min(ndf,nint(pl(1))))
            fnorm = 3
          elseif(jsw.eq.5) then
            isw   = 1
            if(fnorm.eq.3) then ! Error
              write(ilg,3005)
              write(iow,3005)
              call plstop()
            endif
          endif
          do i = 1,3
            errck = pinput(td,5)
            i1 = nint(td(1))
            xl(1,i1) = td(2)
            xl(2,i1) = td(3)
            if(ndm.eq.3) then
              xl(3,i1) = td(4)
            else
              xl(3,i1) = 0.0d0
            endif
            if(isw.ne.2 .and. isw.ne.6) pl(i1)   = td(ndm+2)
          end do ! i

c       Linear segment

        elseif( pcomp(ptype,'line',4)) then
          if(jsw.eq.4) then
            isw   = 1
            ddof  = max(1,min(ndf,nint(pl(1))))
            fnorm = 3
          elseif(jsw.eq.5) then
            isw   = 1
            if(fnorm.eq.3) then ! Error
              write(ilg,3005)
              write(iow,3005)
              call plstop()
            endif
          elseif(jsw.eq.11) then
            call plinlod(mr(np(190)),pl,x,f,prt)
            go to 1
          endif
          do i = 1,2
            errck = pinput(td,5)
            i1 = nint(td(1))
            xl(1,i1) = td(2)
            xl(2,i1) = td(3)
            if(ndm.eq.3) then
              xl(3,i1) = td(4)
            else
              xl(3,i1) = 0.0d0
            endif
            if(isw.ne.2 .and. isw.ne.6) pl(i1)   = td(ndm+2)
          end do ! i
          xl(1,3) = 0.5d0*(xl(1,1) + xl(1,2))
          xl(2,3) = 0.5d0*(xl(2,1) + xl(2,2))
          xl(3,3) = 0.5d0*(xl(3,1) + xl(3,2))
          if(isw.ne.2 .and. isw.ne.6) pl(3) = 0.5d0*(pl(1) + pl(2))
        endif

c       Output data

        if(prt) then
          if(isw.eq.1 .or. isw.eq.5) then
            if(fnorm.le.2) then
              call prtitl(prth)
              write(iow,2000) ((xl(i,j),i=1,2),j=1,3),(pl(i),i=1,3)
              if(ior.lt.0) then
                write(*,2000) ((xl(i,j),i=1,2),j=1,3),(pl(i),i=1,3)
              end if
            elseif(fnorm.eq.3) then
              call prtitl(prth)
              write(iow,2030) ((xl(i,j),i=1,2),j=1,3),(pl(i),i=1,3)
              if(ior.lt.0) then
                write(*,2030) ((xl(i,j),i=1,2),j=1,3),(pl(i),i=1,3)
              end if
            end if
          elseif(isw.eq.2) then
            call prtitl(prth)
            if(ndm.eq.2) then
              write(iow,2010) ((xl(i,j),i=1,2),j=1,3),(i,i=1,ndf)
            elseif(ndm.eq.3) then
              write(iow,2110) ((xl(i,j),i=1,3),j=1,3),(i,i=1,ndf)
            endif
            write(iow,2011) (nint(pl(i)),i=1,ndf)
            if(ior.lt.0) then
              if(ndm.eq.2) then
                write(*,2010) ((xl(i,j),i=1,2),j=1,3),(i,i=1,ndf)
              elseif(ndm.eq.3) then
                write(*,2110) ((xl(i,j),i=1,3),j=1,3),(i,i=1,ndf)
              endif
              write(*,2011) (nint(pl(i)),i=1,ndf)
            end if
          elseif(isw.eq.3) then
            call prtitl(prth)
            write(iow,2020) ((xl(i,j),i=1,2),j=1,3),(pl(i),i=1,3)
            if(ior.lt.0) then
              write(*,2020) ((xl(i,j),i=1,2),j=1,3),(pl(i),i=1,3)
            end if
          elseif(isw.eq.4) then
            call prtitl(prth)
            write(iow,2030) ((xl(i,j),i=1,2),j=1,3),(pl(i),i=1,3)
            if(ior.lt.0) then
              write(*,2030) ((xl(i,j),i=1,2),j=1,3),(pl(i),i=1,3)
            end if
          elseif(isw.eq.6) then
            call prtitl(prth)
            write(iow,2040) ((xl(i,j),i=1,2),j=1,3),(i,i=1,ndf)
            write(iow,2041) (nint(pl(i)),i=1,ndf)
            if(ior.lt.0) then
              write(*,2040) ((xl(i,j),i=1,2),j=1,3),(i,i=1,ndf)
              write(*,2041) (nint(pl(i)),i=1,ndf)
            end if
          elseif(isw.eq.10) then
            call prtitl(prth)
            if(ndm.eq.2) then
              write(iow,2210) ((xl(i,j),i=1,2),j=1,3),(i,i=1,ndf)
            elseif(ndm.eq.3) then
              write(iow,2211) ((xl(i,j),i=1,3),j=1,3),(i,i=1,ndf)
            endif
            write(iow,2011) (nint(pl(i)),i=1,ndf)
            if(ior.lt.0) then
              if(ndm.eq.2) then
                write(*,2210) ((xl(i,j),i=1,2),j=1,3),(i,i=1,ndf)
              elseif(ndm.eq.3) then
                write(*,2211) ((xl(i,j),i=1,3),j=1,3),(i,i=1,ndf)
              endif
              write(*,2011) (nint(pl(i)),i=1,ndf)
            end if
          endif
        endif

c       Transform coordinates for cartesian mode

        if(.not.polfl) then
          do n = 1,3
            xp(1) = xl(1,n)
            xp(2) = xl(2,n)
            xp(3) = xl(3,n)
            xl(1,n) = xr(1) + tr(1,1)*xp(1) + tr(1,2)*xp(2)
     &                      + tr(1,3)*xp(3)
            xl(2,n) = xr(2) + tr(2,1)*xp(1) + tr(2,2)*xp(2)
     &                      + tr(2,3)*xp(3)
            xl(3,n) = xr(3) + tr(3,1)*xp(1) + tr(3,2)*xp(2)
     &                      + tr(3,3)*xp(3)
          end do ! n
        endif

c       Project nodes to points

        errck = .false.
        tol   =  tol0
c       tol   = 3.d0

100     call prj2dl(gap0,tol0,tol,x0,ip,xin,x,xl,ndm,numnp,polfl)

c       Check if end points found

        ximin = +1.d0
        ximax = -1.d0
        do n = 1,numnp
          if(xin(n).ne.0.0d0) then
            ximin = min(ximin,abs(xin(n)))
            ximax = max(ximax,abs(xin(n)))
          end if
        end do ! n

c       For cylindrical coordinates check number of 1 and -1 values

        if(polfl) then
          i1 = 0
          i2 = 0
          do n = 1,numnp
            if(xin(n).ge.  1.d0 - tol0) then
              i1 = i1 + 1
            endif
            if(xin(n).le. -1.d0 + tol0) then
              i2 = i2 + 1
            endif
          end do ! n
          if(i1.ge.2) then
            ximin = -1.d0
          endif
          if(i2.ge.2) then
            ximax =  1.d0
          endif
        endif

        tol = tol0 + max(1.d0+ximin,1.d0-ximax,0.0d0)
c       tol = tol0 + max(1.d0-ximin,0.0d0)

        errck = .not.errck
        errck = .false.
        if(errck .and. tol.gt.tol0+tolxi) go to 100

c       Loop through elements to find adjacent nodes

        do n = 1,numel
          ep(n,1) = 0
          ep(n,2) = 0
          ep(n,3) = 0
          ep(n,4) = 0

c         Check number of nodes on element

          do j = 1,nen
            if(ix(j,n).ne.0) then
              nel = j
            endif
          end do ! j

          if(ndm.eq.2 .and. (nel.eq.6 .or. nel.eq.7)) then
            jmin = 3
          else
            jmin = min(4,nel)
          endif

          do j = 1,jmin
            i1 = ix(j,n)
            if(i1.ne.0) then
              if(ip(i1).ne.0) then
                if(j.lt.jmin) then
                  i2 = ix(j+1,n)
                  if(i2.eq.0) then
                    i2 = ix(1,n)
                  endif
                else
                  i2 = ix(1,n)
                endif
                if(ip(i2).gt.0 .and. i1.ne.i2) then
c                 if(ep(n,1).le.0 .and. ep(n,2).le.0 ) then
                  if(ep(n,1).le.0 .and. ep(n,2).le.0 .and.
     &              (min(abs(xin(i1)),abs(xin(i2))).le.1.d0+tolxi)) then
                    ep(n,1) = i1
                    ep(n,2) = i2

c                   Triangles

                    if(nel.eq.6 .or. nel.eq.7) then

                      i3 = ix(j+3,n)
                      if(i3.gt.0) then
                        if(ip(i3).gt.0) ep(n,3) = i3
                      endif

c                   Cubic elements

                    elseif(nel.eq.16) then

                      i3 = ix(2*j+3,n)
                      if(i3.gt.0) then
                        if(ip(i3).gt.0) ep(n,3) = i3
                      endif
                      i3 = ix(2*j+4,n)
                      if(i3.gt.0) then
                        if(ip(i3).gt.0) ep(n,4) = i3
                      endif

c                   Quadratic elements

                    elseif(j+4.le.nel) then

                      i3 = ix(j+4,n)
                      if(i3.gt.0) then
                        if(ip(i3).gt.0) ep(n,3) = i3
                      endif
                      ep(n,4) = 0
                    endif

                  endif
                endif
              endif
            endif
          end do ! j
        end do ! n

c       Remove duplicates

        do n = 1,numel-1
          if(ep(n,1).gt.0) then
            do m = n+1,numel
              if(ep(m,1).eq.ep(n,1) .and. ep(m,2).eq.ep(n,2) .or.
     &           ep(m,1).eq.ep(n,2) .and. ep(m,2).eq.ep(n,1)) then
                ep(m,1) = 0
                ep(m,2) = 0
                ep(m,3) = 0
                ep(m,4) = 0
              end if
            end do ! m
          end if
        end do ! n

c       Check for polar projections: Possible that angle 0 and 360 both
c       project to 'xi' of -1.0d0

        if(isw.eq.1 .and. polfl) then
          do n = 1,numel
            if(ep(n,1).gt.0) then
              xia = xin(ep(n,1))
              xib = xin(ep(n,2))
              if(xib.lt.xia) then
                if(abs(xib).ge.1.d0-tol0) then
                  xin(ep(n,2)) = abs(xib)
                else
                  write(ilg,3002) n,ep(n,1),ep(n,2)
                  write(iow,3002) n,ep(n,1),ep(n,2)
                  if(ior.lt.0) then
                    write(*,3002) n,ep(n,1),ep(n,2)
                  endif
                endif
              endif
            endif
          end do ! n
        endif

c       Compute nodal forces from surface pressures

        if(isw.eq.1 .and. fnorm.le.2) then

c         Find end points

          do n = 1,numnp
            ip(n) = 0
          end do ! n

          xld(1) = 0.5d0*(xl(1,2) - xl(1,1))
          xld(2) = 0.5d0*(xl(2,2) - xl(2,1))
          xqd(1) = xl(1,1) + xl(1,2) - 2.d0*xl(1,3)
          xqd(2) = xl(2,1) + xl(2,2) - 2.d0*xl(2,3)
          do n = 1,numel
            if(ep(n,1).gt.0) then

c             Order segments along master direction

              if(polfl) then
                yc = xl(2,3) + xin(ep(n,1))*(xld(2)
     &                       + 0.5d0*xin(ep(n,1))*xqd(2))
                if(xl(2,2).gt.xl(2,1)) then
                  call pdegree(yc, cn,sn)
                  cn = -cn
                else
                  call pdegree(yc, cn,sn)
                  sn = -sn
                endif
              else
                cn = xld(1) + xin(ep(n,1))*xqd(1)
                sn = xld(2) + xin(ep(n,1))*xqd(2)
              endif

              d = (x(1,ep(n,2)) - x(1,ep(n,1)))*cn
     &          + (x(2,ep(n,2)) - x(2,ep(n,1)))*sn

              if(d.lt.0.0d0) then
                ep(n,1) = 0
                ep(n,2) = 0
                ep(n,3) = 0
                ep(n,4) = 0
              end if

c             Count occurances of node

              if(ep(n,1).gt.0) ip(ep(n,1)) = ip(ep(n,1)) + 1
              if(ep(n,2).gt.0) ip(ep(n,2)) = ip(ep(n,2)) + 1
            end if
          end do ! n

c         Set end point array

          j         = 0
          nend(1,1) = 0
          nend(2,1) = 0
          do n = 1,numnp
            if(ip(n).eq.1) then
              j         = j + 1
              nend(j,1) = n
            end if
          end do ! n

c         Check for error

          if(nend(1,1).eq.0 .or. nend(2,1).eq.0) then
            write(iow,3004)
            write(ilg,3004)
            if(ior.lt.0) then
              write(*,3004)
            endif
          endif

c         Compute nodal forces for pressures

          if(prt) then
            write(iow,2001)
            if(ior.lt.0) write(*,2001)
          endif

          do n = 1,numel

            i1 = ep(n,1)

            if(i1.gt.0) then

              i2     = ep(n,2)
              xia    = xin(i1)
              xib    = xin(i2)
              if(ep(n,3).eq.0) then
                dxl(1) = x(1,i2) - x(1,i1)
                dxl(2) = x(2,i2) - x(2,i1)

                if(xia.lt.-1.d0) then
                  if(xib.gt.1.0d0) then
                    pa    = pl(1)
                    pb    = pl(2)
                    pc    = pl(3)
                    alpha = 1.d0/(xib - xia)
                    pa    = alpha*(pa + pc + pc)*one3
                    pb    = alpha*(pc + pc + pb)*one3
                    fb    = ((-1.d0 - xia)*pa + (1 - xia)*pb)*alpha
                    fa    = pa + pb - fb
                  elseif (xib.gt.-1.d0) then
                    xic   = -0.5d0 + xib*0.5d0
                    pa    = pl(1)
                    pb    = pl(1)*0.5d0*xib*(xib-1.d0)
     &                    + pl(3)*(1.d0-xib*xib)
     &                    + pl(2)*0.5d0*xib*(xib+1.d0)
                    pc    = pl(1)*0.5d0*xic*(xic-1.d0)
     &                    + pl(3)*(1.d0-xic*xic)
     &                    + pl(2)*0.5d0*xic*(xic+1.d0)
                    if(xib-xia .gt. tolxi) then
                      alpha = (xib + 1.0d0)/(xib - xia)
                    else
                      alpha = 1.d0
                    endif
                    fa    = (pa + pc + pc)*one6
                    fb    = (pc + pc + pb)*one6
                    fb    = alpha*(fa + fb)
                    fa    = alpha*alpha*fa
                    fb    = fb - fa
                  else
                    pa    = 0.0d0
                    pb    = 0.0d0
                    pc    = 0.0d0
                    fa    = 0.0d0
                    fb    = 0.0d0
                  end if
                elseif (xia.lt.1.d0) then

                  pa    = pl(1)*0.5d0*xia*(xia-1.d0)
     &                  + pl(3)*(1.d0-xia*xia)
     &                  + pl(2)*0.5d0*xia*(xia+1.d0)

                  if(xib.lt.1.d0) then

                    xic   = (xia + xib)*0.5d0
                    pb    = pl(1)*0.5d0*xib*(xib-1.d0)
     &                    + pl(3)*(1.d0-xib*xib)
     &                    + pl(2)*0.5d0*xib*(xib+1.d0)
                    alpha = 1.d0

                  else

                    xic   = xia*0.5d0 + 0.5d0
                    pb    = pl(2)
                    if(xib-xia .gt. tolxi) then
                      alpha = (1.d0 - xia)/(xib - xia)
                    else
                      alpha = 1.d0
                    endif

                  end if
                  pc    = pl(1)*0.5d0*xic*(xic-1.d0)
     &                  + pl(3)*(1.d0-xic*xic)
     &                  + pl(2)*0.5d0*xic*(xic+1.d0)

                  fa    = (pa + pc + pc)*one6
                  fb    = (pc + pc + pb)*one6
                  fa    = alpha*(fa + fb)
                  fb    = alpha*alpha*fb
                  fa    = fa - fb
                else
                  fa = 0.0d0
                  fb = 0.0d0
                  pa = 0.0d0
                  pb = 0.0d0
                end if

c               Add forces

                if(abs(fa)+abs(fb).gt.0.0d0) then

                  if(fnorm.eq.1) then
                    fxa =  fa*dxl(2)
                    fya = -fa*dxl(1)
                    fxb =  fb*dxl(2)
                    fyb = -fb*dxl(1)
                  elseif(fnorm.eq.2) then
                    fxa =  fa*dxl(1)
                    fya =  fa*dxl(2)
                    fxb =  fb*dxl(1)
                    fyb =  fb*dxl(2)
                  endif

c                 Check for sloping boundary conditions

                  if(ang(i1).ne.0.0d0) then
                    call pdegree(ang(i1), sn,cn)
                    ff  =  cn*fxa + sn*fya
                    fya = -sn*fxa + cn*fya
                    fxa =  ff
                  endif

                  if(ang(i2).ne.0.0d0) then
                    call pdegree(ang(i2), sn,cn)
                    ff  =  cn*fxb + sn*fyb
                    fyb = -sn*fxb + cn*fyb
                    fxb =  ff
                  endif

c                 Normal element

                  if(ndf.lt.8) then

                    f(1,i1,1) =  f(1,i1,1) + fxa
                    f(2,i1,1) =  f(2,i1,1) + fya
                    f(1,i2,1) =  f(1,i2,1) + fxb
                    f(2,i2,1) =  f(2,i2,1) + fyb

c                 Hierarchic formulation

                  else

                    lint = 3
                    call int1d(lint,sw)
                    do l = 1,lint
                      call shap1d(sw(1,l),2,shp)
                      dxl(1) = shp(1,1)*x(1,i1) + shp(1,2)*x(1,i2)
                      dxl(2) = shp(1,1)*x(2,i1) + shp(1,2)*x(2,i2)
                      xc     = shp(2,1)*x(1,i1) + shp(2,2)*x(1,i2)
                      yc     = shp(2,1)*x(2,i1) + shp(2,2)*x(2,i2)
                      ff = (pa*0.5d0*sw(1,l)*(sw(1,l)-1.d0)
     &                   +  pc*(1.d0-sw(1,l)*sw(1,l))
     &                   +  pb*0.5d0*sw(1,l)*(sw(1,l)+1.d0))*sw(2,l)
                      x1        = xc - x(1,i1)
                      y1        = yc - x(2,i1)
                      x2        = x1*x1
                      xy        = x1*y1
                      y2        = y1*y1
                      f(1,i1,1) = f(1,i1,1) + dxl(2)*shp(2,1)*ff
                      f(2,i1,1) = f(2,i1,1) - dxl(1)*shp(2,1)*ff
                      f(3,i1,1) = f(3,i1,1) + x2*dxl(2)*shp(2,1)*ff
                      f(4,i1,1) = f(4,i1,1) - x2*dxl(1)*shp(2,1)*ff
                      f(5,i1,1) = f(5,i1,1) + xy*dxl(2)*shp(2,1)*ff
                      f(6,i1,1) = f(6,i1,1) - xy*dxl(1)*shp(2,1)*ff
                      f(7,i1,1) = f(7,i1,1) + y2*dxl(2)*shp(2,1)*ff
                      f(8,i1,1) = f(8,i1,1) - y2*dxl(1)*shp(2,1)*ff

                      x1        = xc - x(1,i2)
                      y1        = yc - x(2,i2)
                      x2        = x1*x1
                      xy        = x1*y1
                      y2        = y1*y1
                      f(1,i2,1) = f(1,i2,1) + dxl(2)*shp(2,2)*ff
                      f(2,i2,1) = f(2,i2,1) - dxl(1)*shp(2,2)*ff
                      f(3,i2,1) = f(3,i2,1) + x2*dxl(2)*shp(2,2)*ff
                      f(4,i2,1) = f(4,i2,1) - x2*dxl(1)*shp(2,2)*ff
                      f(5,i2,1) = f(5,i2,1) + xy*dxl(2)*shp(2,2)*ff
                      f(6,i2,1) = f(6,i2,1) - xy*dxl(1)*shp(2,2)*ff
                      f(7,i2,1) = f(7,i2,1) + y2*dxl(2)*shp(2,2)*ff
                      f(8,i2,1) = f(8,i2,1) - y2*dxl(1)*shp(2,2)*ff
                    end do ! l

                  endif

                  if(prt) then
                    write(iow,2002) i1,fxa,fya,i2,fxb,fyb
                    if(ior.lt.0) then
                      write(*,2002) i1,fxa,fya,i2,fxb,fyb
                    endif
                  endif

                end if

c             Numerically integrated quadratic and cubic edge loading

              else

c               Set order

                if(ep(n,4).ne.0) then
                  nel = 4
                else
                  nel = 3
                endif

c               Compute end values for each coordinate

                do i = 1,nel
                  xe(1,i) = x(1,ep(n,i))
                  xe(2,i) = x(2,ep(n,i))
                  fl(1,i) = 0.0d0
                  fl(2,i) = 0.0d0
                end do ! i

                zeta1 = xia
                zeta2 = xib

                xi1 = max(-1.d0, -(2.d0+zeta1+zeta2)/(zeta2-zeta1))
                xi2 = min( 1.d0,  (2.d0-zeta1-zeta2)/(zeta2-zeta1))

                df  = 0.5d0*(xi2 - xi1)

                lint = nel
                call int1d(lint,sw)
                do l = 1,lint

c                 Coordinates for facet and interpolation surface

                  xi   = 0.5d0*((1.d0 - sw(1,l))*xi1
     &                        + (1.d0 + sw(1,l))*xi2)
                  zeta = 0.5d0*((1.d0 - xi)*zeta1  + (1.d0 + xi)*zeta2)

c                 Magnitude of loading

                  ff    = (pl(1)*0.5d0*zeta*(zeta-1.d0)
     &                  +  pl(2)*0.5d0*zeta*(zeta+1.d0)
     &                  +  pl(3)*(1.d0-zeta*zeta))*sw(2,l)*df

c                 Shape functions for facet

                  call shap1d(xi,nel,shp)

c                 Tangent vector to facet

                  xp(1) = shp(1,1)*xe(1,1)
                  xp(2) = shp(1,1)*xe(2,1)
                  do i = 2,nel
                    xp(1) = xp(1) + shp(1,i)*xe(1,i)
                    xp(2) = xp(2) + shp(1,i)*xe(2,i)
                  end do ! i

c                 Tangential loading values

                  do i = 1,nel
                    fl(1,i) = fl(1,i) + shp(2,i)*xp(1)*ff
                    fl(2,i) = fl(2,i) + shp(2,i)*xp(2)*ff
                  end do ! i
                end do ! l

c               Transform for normal loading

                if(fnorm.eq.1) then
                  do i = 1,nel
                    ff      =  fl(1,i)
                    fl(1,i) =  fl(2,i)
                    fl(2,i) = -ff
                  end do ! i
                endif

c               Perform outputs and add to nodal forces

                do i = 1,nel
                  i1 = ep(n,i)

                  if(prt) then
                    write(iow,2002) i1,fl(1,i),fl(2,i)
                    if(ior.lt.0) then
                      write(*,2002) i1,fl(1,i),fl(2,i)
                    endif
                  endif

c                 Sloping boundary conditions

                  if(ang(i1).ne.0.0d0) then
                    call pdegree(ang(i1), sn,cn)
                    ff      =  cn*fl(1,i) + sn*fl(2,i)
                    fl(2,i) = -sn*fl(1,i) + cn*fl(2,i)
                    fl(1,i) =  ff
                  endif

c                 Add to nodal forces

                  f(1,i1,1) =  f(1,i1,1) + fl(1,i)
                  f(2,i1,1) =  f(2,i1,1) + fl(2,i)
                end do ! i

              end if

            end if

          end do ! n

c       Other options

        else

          do n = 1,numnp
            ip(n) = 0
          end do ! n

          do n = 1,numel
            if(ep(n,1).gt.0) then
              ip(ep(n,1)) = 1
              ip(ep(n,2)) = 1
              if(ep(n,3).gt.0) then
                ip(ep(n,3)) = 1
              endif
              if(ep(n,4).gt.0) then
                ip(ep(n,4)) = 1
              endif
            end if
          end do ! n

c         Set specified displacements on line

          if(isw.eq.1 .and. fnorm.eq.3) then

            if(prt) then
              write(iow,2031)
              if(ior.lt.0) write(*,2031)
            end if

            do n = 1,numnp
              if(ip(n).gt.0 .and. abs(xin(n)).le.1.d0+tolxi) then
                xia    = xin(n)
                f(ddof,n,2) = pl(1)*0.5d0*xia*(xia-1.d0)
     &                      + pl(3)*(1.d0-xia*xia)
     &                      + pl(2)*0.5d0*xia*(xia+1.d0)
                if(prt) then
                  write(iow,2032) n,ddof,f(ddof,n,2)
                  if(ior.lt.0) then
                    write(*,2032) n,ddof,f(ddof,n,2)
                  endif
                endif

              endif
            end do ! n

c         Set boundary conditions for nodes on surface

          elseif(isw.eq.2) then

            if(prt) then
              write(iow,2012) (i,i=1,ndf)
              if(ior.lt.0) write(*,2012) (i,i=1,ndf)
            end if
            do n = 1,numnp
              if(ip(n).gt.0 .and. abs(xin(n)).le.1.d0+tolxi) then
                do i = 1,ndf
                  if(nint(pl(i)).ge.0) then
                    id(i,n,2) = id(i,n,2) + nint(pl(i))
                  else
                    id(i,n,2) = 0
                  endif
                end do ! i
                if(prt) then
                  write(iow,2013) n,(id(i,n,2),i=1,ndf)
                  if(ior.lt.0) then
                    write(*,2013) n,(id(i,n,2),i=1,ndf)
                  endif
                endif
              endif
            end do ! n

c         Set boundary angles for nodes on surface

          elseif(isw.eq.3) then

            if(prt) then
              write(iow,2021)
              if(ior.lt.0) write(*,2021)
            endif

            do n = 1,numnp
              if(ip(n).gt.0) then
                xia    = xin(n)
                ang(n) = pl(1)*0.5d0*xia*(xia-1.d0)
     &                 + pl(3)*(1.d0-xia*xia)
     &                 + pl(2)*0.5d0*xia*(xia+1.d0)
                if(prt) then
                  write(iow,2002) n,ang(n)
                  if(ior.lt.0) then
                    write(*,2002) n,ang(n)
                  endif
                endif

              endif
            end do ! n

c         Set displacements for nodes on surface

          elseif(isw.eq.4) then

            if(prt) then
              write(iow,2042) (i,i=1,ndf)
              if(ior.lt.0) then
                write(*,2042) (i,i=1,ndf)
              endif
            endif

            do n = 1,numnp
              if(ip(n).gt.0) then
                do i = 1,ndf
                  f(i,n,2) = pl(i)
                end do ! i
                if(prt) then
                  write(iow,2043) n,(f(i,n,2),i=1,ndf)
                  if(ior.lt.0) then
                    write(*,2043) n,(f(i,n,2),i=1,ndf)
                  endif
                endif

              endif

            end do ! n

c         Set proportional load numbers for surface

          elseif(isw.eq.6) then

            if(prt) then
              write(iow,2052)
              if(ior.lt.0) write(*,2052)
            endif

            do n = 1,numnp
              if(ip(n).gt.0) then
                fp(1) = np(29) + ndf*(n-1) - 1
                do i = 1,ndf
                  mr(fp(1)+i) = nint(pl(i))
                end do ! i
                if(prt) then
                  write(iow,2053) n,(mr(fp(1)+i),i=1,ndf)
                  if(ior.lt.0) then
                    write(*,2053) n,(mr(fp(1)+i),i=1,ndf)
                  endif
                endif

              endif
            end do ! n

c         Set base conditions for nodes on surface

          elseif(isw.eq.10) then

            if(prt) then
              write(iow,2212) (i,i=1,ndf)
              if(ior.lt.0) write(*,2212) (i,i=1,ndf)
            end if
            do n = 1,numnp
              fp(1) = np(125) + ndf*(n-1) - 1
              if(ip(n).gt.0 .and. abs(xin(n)).le.1.d0+tolxi) then
                do i = 1,ndf
                  if(nint(pl(i)).ge.0) then
                    mr(fp(1)+i) = mr(fp(1)+i) + nint(pl(i))
                  else
                    mr(fp(1)+i) = 0
                  endif
                end do ! i
                if(prt) then
                  write(iow,2013) n,(mr(fp(1)+i),i=1,ndf)
                  if(ior.lt.0) then
                    write(*,2013) n,(mr(fp(1)+i),i=1,ndf)
                  endif
                endif
              endif
            end do ! n

          endif

        endif

        go to 1

c       End of current file

300     continue

      end do ! nn

      prt  = oprt
      prth = oprth

c     Formats

2000  format('   C o o r d i n a t e    S u r f a c e   L o a d s'/
     &       '   x_1 coord.   y_1 coord.   x_2 coord.   y_2 coord.',
     &       '   x_3 coord.   y_3 coord.'/1p,6e13.4/
     &       '   p_1 press.                p_2 press.             ',
     &       '   p_3 press.'/3(1p,1e13.4,13x))

2001  format(/'       N o d a l    F o r c e s'//
     &       '       a_node    x_a force    y_a force       b_node',
     &       '    x_b force    y_b force'/)

2002  format(2(i13,1p,2e13.4))

2004  format(/'       P o l a r    S u r f a c e'//
     &       '     x_coord.     y_coord.     z_coord.'/1p,3e13.4)

2010  format('   C o o r d i n a t e    S u r f a c e',
     &       '   C o n d i t i o n'/
     &       '   x_1 coord.   y_1 coord.   x_2 coord.   y_2 coord.',
     &       '   x_3 coord.   y_3 coord.'/1p,6e13.4/
     &       (10(i3,' B.C.',:)))

2110  format('   C o o r d i n a t e    S u r f a c e',
     &       '   C o n d i t i o n'/
     &       '   x_1 =',1p,1e13.4,'   y_1 =',1p,1e13.4,
     &       '   z_1 =',1p,1e13.4/
     &       '   x_2 =',1p,1e13.4,'   y_2 =',1p,1e13.4,
     &       '   z_2 =',1p,1e13.4/
     &       '   x_3 =',1p,1e13.4,'   y_3 =',1p,1e13.4,
     &       '   z_3 =',1p,1e13.4/
     &       (10(i3,' B.C.',:)))

2011  format(10i8)

2012  format(/'       N o d a l    B o u n d a r y    C o d e s'//
     &       '      Node',6(i5,'-B.C.':)/(10x,6(i5,'-B.C.':)))

2013  format(7i10/(10x,6i10))

2020  format('   C o o r d i n a t e    S u r f a c e   L o a d s'/
     &       '   x_1 coord.   y_1 coord.   x_2 coord.   y_2 coord.',
     &       '   x_3 coord.   y_3 coord.'/1p,6e13.4/
     &       '   a_1 angle                 a_2 angle              ',
     &       '   a_3 angle '/3(1p,1e13.4,13x))

2021  format(/'       N o d a l    A n g l e s'//
     &       '       a_node    t_a angle')

2030  format('   C o o r d i n a t e    S u r f a c e   D i s p l.'/
     &       '   x_1 coord.   y_1 coord.   x_2 coord.   y_2 coord.',
     &       '   x_3 coord.   y_3 coord.'/1p,6e13.4/
     &       '   d_1 displ                 d_2 displ              ',
     &       '   d_3 displ '/3(1p,1e13.4,13x))

2031  format(/'       N o d a l    D i s p l a c e m e n t s'//
     &       '       a_node    a_dof   d_a displ')

2032  format(i13,i8,1p,e13.4)

2040  format('   C o o r d i n a t e    S u r f a c e',
     &       '   C o n d i t i o n'/
     &       '   x_1 coord.   y_1 coord.   x_2 coord.   y_2 coord.',
     &       '   x_3 coord.   y_3 coord.'/1p,6e13.4/
     &       (10(i3,' Prop',:)))

2041  format(10i8)

2042  format(/'       N o d a l    B o u n d a r y    D i s p l.'//
     &       '    Node',6(i6,'-Displ':)/(10x,6(i6,'-Displ':)))

2043  format(1i8,1p,6e12.4/(8x,1p,6e12.4))

2052  format(/'       N o d a l    P r o p o r t i o n a l    ',
     &        'L o a d s'//
     &        '      Node',6(i3,'-Propld':)/(10x,6(i3,'-Propld':)))

2053  format(7i10/(10x,6i10))

2210  format('   C o o r d i n a t e    S u r f a c e',
     &       '   C o n d i t i o n'/
     &       '   x_1 coord.   y_1 coord.   x_2 coord.   y_2 coord.',
     &       '   x_3 coord.   y_3 coord.'/1p,6e13.4/
     &       (10(i3,' Base',:)))

2211  format('   C o o r d i n a t e    S u r f a c e',
     &       '   C o n d i t i o n'/
     &       '   x_1 =',1p,1e13.4,'   y_1 =',1p,1e13.4,
     &       '   z_1 =',1p,1e13.4/
     &       '   x_2 =',1p,1e13.4,'   y_2 =',1p,1e13.4,
     &       '   z_2 =',1p,1e13.4/
     &       '   x_3 =',1p,1e13.4,'   y_3 =',1p,1e13.4,
     &       '   z_3 =',1p,1e13.4/
     &       (10(i3,' Base',:)))

2212  format(/'       N o d a l    B a s e    C o d e s'//
     &       '      Node',6(i5,'-base':)/(10x,6(i5,'-base':)))


3001  format(' *ERROR* PESURF: Segment node too large, node =',i5)

3002  format(' *ERROR* PESURF: Segment for element',i8,
     &       ' has bad projection values for nodes',i8,' and',i8)

3003  format(' *ERROR* PESURF: In surface condition inputs')

3004  format(' *ERROR* PESURF: No surface located')

3005  format(' *ERROR* PESURF: Attempt to specify DISPlacement during',
     &       ' CFORce input.')

      end
