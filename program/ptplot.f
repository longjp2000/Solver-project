c$Id: ptplot.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine ptplot (x,ix,ct,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set tplot options

c      Inputs:
c         x(ndm,*)   - Nodal coordinates
c         ix(nen1,*) - Element connection list
c         ct(3)      - Command parameters for current command
c         prt        - Flag, print data if true

c      Outputs:
c         Data lists for tplot
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'codat.h'
      include  'endata.h'
      include  'fdata.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'ioincl.h'
      include  'iosave.h'
      include  'part0.h'
      include  'part3.h'
      include  'part7.h'
      include  'pbody.h'
      include  'plcapt.h'
      include  'plist.h'
      include  'prflag.h'
      include  'prlod.h'
      include  'ptdat1.h'
      include  'ptdat2.h'
      include  'ptdat3.h'
      include  'ptdat4.h'
      include  'ptdat5.h'
      include  'ptdat6.h'
      include  'ptdat7.h'
      include  'ptdat8.h'
      include  'ptdat9.h'
      include  'ptdata.h'
      include  'ptdatb.h'
      include  'sdata.h'
      include  'tdata.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   prt, pcomp, errck, tinput
      character type*4,yyy*15
      integer   i, j, k, ii, jj, kk, n, nn
      integer   npmx,nsmx
      real*8    dotx, dd, dist, ct(3),td(5),xc(3)
      integer   ix(nen1,*)
      real*8    x(ndm,*)

      save

      data      npmx /200/, nsmx/200/

c     [tplo]t,,<interval> - time history plots for displacements / stresses

c      Options : disp n1 n2 x y z
c                velo n1 n2 x y z
c                acce n1 n2 x y z
c                reac n1 n2 x y z
c                stre n1 n2 x y z
c                elem n1 n2 x y z
c                user n1 n2 x y z
c                cont n1 n2 x y z
c                arcl n1 n2
c                rsum n1 n2
c                sums n1 n2
c                ener
c                show

c     Set output incrment

      ntincr = max(1,int(ct(1)))

      if(prt) then
        write(iow,2001) ntincr
        if(ior.lt.0) then
          write(*,2001) ntincr
        endif
      endif

c     Input list of time history output quantities

      if(ior.lt.0) write(*,3001)
21    if(ior.lt.0) call pprint('          >')
      errck = tinput(yyy,1,td,5)
      if(errck) go to 21
      type = yyy
      n    = td(1)
      i    = td(2)

c     Find the location for unspecified node/elements

      if(n.eq.0) then

c       Locate a node

        if(pcomp(yyy,'disp',4) .or.
     &     pcomp(yyy,'velo',4) .or.
     &     pcomp(yyy,'acce',4) .or.
     &     pcomp(yyy,'reac',4) .or.
     &     pcomp(yyy,'cont',4) .or.
     &     pcomp(yyy,'arcl',4)) then

          dist = 0.d+0
          do nn = 1,numnp
            do n = 1,ndm
              dist = max(dist,abs(x(n,nn)))
            end do ! n
          end do ! n

          do nn = 1,numnp
            if(mr(np(190)+nn-1).ge.0) then
              dd = dotx(x(1,nn),td(3),ndm)
              dd = sqrt(dd)
              if(dd.lt.dist) then
                dist = dd
                n    = nn
              endif
            endif
          end do ! nn

c       Locate an element

        elseif(pcomp(yyy,'stre',4) .or. pcomp(yyy,'elem',4)) then

          dist = 0.d+0
          do nn = 1,numnp
            do n = 1,ndm
              dist = max(dist,abs(x(n,nn)))
            end do ! n
          end do ! n

          do nn = 1,numel
            do jj = 1,ndm
              xc(jj) = 0.0d0
            end do ! jj

            kk = 0
            do jj = 1,nen
              k = ix(jj,nn)
              if(k.gt.0) then
                kk    = kk + 1
                do ii = 1,ndm
                  xc(ii) = xc(ii) + x(ii,k)
                end do ! ii
              endif
            end do ! jj
            if(kk.gt.0) then
              do ii = 1,ndm
                xc(ii) = xc(ii) / kk
              end do ! ii

              dd = dotx(xc,td(3),ndm)
              dd = sqrt(dd)
              if(dd.lt.dist) then
                dist = dd
                n    = nn
              endif
            endif
          end do ! nn

        endif
      endif

      if(pcomp(type,'    ',4)) return

c     Displacements

      if(pcomp(type,'disp',4)) then
        ndplts         = min(npmx,ndplts + 1)
        idpl(1,ndplts) = ndf*(n-1)+i
        idpl(2,ndplts) = n

c     Velocities

      elseif(pcomp(type,'velo',4)) then
        nvplts         = min(npmx,nvplts + 1)
        ivpl(1,nvplts) = ndf*(n-1)+i
        ivpl(2,nvplts) = n

c     Accelerations

      elseif(pcomp(type,'acce',4)) then
        naplts         = min(npmx,naplts + 1)
        iapl(1,naplts) = ndf*(n-1)+i
        iapl(2,naplts) = n

c     Stresses

      elseif(pcomp(type,'stre',4) .or. pcomp(type,'elem',4)) then
        nsplts         = min(nsmx,nsplts + 1)
        ispl(1,nsplts) = n
        ispl(2,nsplts) = i

c     User Stresses

      elseif(pcomp(type,'user',4)) then
        nuplts         = min(npmx,nuplts + 1)
        iupl(1,nuplts) = n
        iupl(2,nuplts) = i

c     Reactions

      elseif(pcomp(type,'reac',4)) then
        nrplts         = min(npmx,nrplts + 1)
        irpl(1,nrplts) = ndf*(n-1)+i
        irpl(2,nrplts) = n

c     Energy

      elseif(pcomp(type,'ener',4)) then
        neplts         = min(npmx,neplts + 1)
        iepl(1,neplts) = n
        iepl(2,neplts) = i

c     Contacts

      elseif(pcomp(type,'cont',4)) then
        ncplts         = min(npmx,ncplts + 1)
        icpl(1,ncplts) = n
        icpl(2,ncplts) = i

c     Arclength

      elseif(pcomp(type,'arcl',4)) then
        nlplts         = min(npmx,nlplts + 1)
        ilpl(1,nlplts) = ndf*(n-1)+i
        ilpl(2,nlplts) = n

c     Reaction sum

      elseif(pcomp(type,'rsum',4)) then
        nqplts         = min(npmx,nqplts + 1)
        iqpl(1,nqplts) = n            ! Dof to sum
        iqpl(2,nqplts) = i            ! Initial node
        iqpl(3,nqplts) = nint(td(3))  ! Final   node

c     Sums

      elseif(pcomp(type,'sums',4)) then
        ntplts         = min(npmx,ntplts + 1)
        itpl(1,ntplts) = n
        itpl(2,ntplts) = i
        tpld(1,ntplts) = td(3)
        tpld(2,ntplts) = td(4)
        tpld(3,ntplts) = td(5)

c     Material states

      elseif(pcomp(type,'mate',4)) then
        nmplts         = min(10,nmplts + 1)
        impl(1,nmplts) = max(1,n)
        impl(2,nmplts) = i

c     Show: Active outputs

      elseif(pcomp(type,'show',4)) then

        do n = 1,ndplts
          k = idpl(2,n)
          i = idpl(1,n) - ndf*(k -1)
          if(ior.lt.0) then
            write(*,3003) n,idpl(2,n),i,(x(j,k),j=1,ndm)
          endif
          write(iow,3003) n,idpl(2,n),i,(x(j,k),j=1,ndm)
        end do ! n

        do n = 1,nvplts
          k = ivpl(2,n)
          i = ivpl(1,n) - ndf*(k -1)
          if(ior.lt.0) then
            write(*,3004) n,ivpl(2,n),i,(x(j,k),j=1,ndm)
          endif
          write(iow,3004) n,ivpl(2,n),i,(x(j,k),j=1,ndm)
        end do ! n

        do n = 1,naplts
          k = iapl(2,n)
          i = iapl(1,n) - ndf*(k -1)
          if(ior.lt.0) then
            write(*,3005) n,iapl(2,n),i,(x(j,k),j=1,ndm)
          endif
          write(iow,3005) n,iapl(2,n),i,(x(j,k),j=1,ndm)
        end do ! n

        do n = 1,nsplts
          if(ior.lt.0) then
            write(*,3006) n,ispl(1,n),ispl(2,n)
          endif
          write(iow,3006) n,ispl(1,n),ispl(2,n)
        end do ! n

        do n = 1,ncplts
          if(ior.lt.0) then
            write(*,3007) n,icpl(1,n),icpl(2,n)
          endif
          write(iow,3007) n,icpl(1,n),icpl(2,n)
        end do ! n

        do n = 1,nrplts
          k = irpl(2,n)
          i = irpl(1,n) - ndf*(k -1)
          if(ior.lt.0) then
            write(*,3008) n,irpl(2,n),i,(x(j,k),j=1,ndm)
          endif
          write(iow,3008) n,irpl(2,n),i,(x(j,k),j=1,ndm)
        end do ! n

        do n = 1,neplts
          if(ior.lt.0) then
            write(*,3009) n,iepl(1,n),iepl(2,n)
          endif
          write(iow,3009) n,iepl(1,n),iepl(2,n)
        end do ! n

        do n = 1,nlplts
          i = ilpl(1,n) - ndf*(ilpl(2,n) -1)
          if(ior.lt.0) then
            write(*,3010) n,ilpl(2,n),i
          endif
          write(iow,3010) n,ilpl(2,n),i
        end do ! n

        do n = 1,ntplts
          if(ior.lt.0) then
            write(*,3011) n,itpl(1,n),itpl(2,n),tpld(1,n),tpld(2,n)
          endif
          write(iow,3011) n,itpl(1,n),itpl(2,n),tpld(1,n),tpld(2,n)
        end do ! n

        do n = 1,nmplts
          if(ior.lt.0) then
            write(*,3012) n,impl(1,n),impl(2,n)
          endif
          write(iow,3012) n,impl(1,n),impl(2,n)
        end do ! n

        do n = 1,nuplts
          if(ior.lt.0) then
            write(*,3013) n,iupl(1,n),iupl(2,n)
          endif
          write(iow,3013) n,iupl(1,n),iupl(2,n)
        end do ! n

        do n = 1,nqplts
          if(ior.lt.0) then
            write(*,3014) n,iqpl(1,n),iqpl(2,n),iqpl(3,n)
          endif
          write(iow,3014) n,iqpl(1,n),iqpl(2,n),iqpl(3,n)
        end do ! n

        return

      endif
      go to 21

c     Formats

2001  format(/'   Output interval for time history data =',i4)

3001  format(' Input: Type (disp:velo:acce:stres:cont:reac:arcl:ener:',
     &       'sums,user);'/'        Node/Elmt; dof/no.')

3003  format(1x,'Plot',i3,' Displ. : Node  =',i8,' DOF =',i3,
     &       ' X =',1p,3e11.3)

3004  format(1x,'Plot',i3,' Veloc. : Node  =',i8,' DOF =',i3,
     &       ' X =',1p,3e11.3)

3005  format(1x,'Plot',i3,' Accel. : Node  =',i8,' DOF =',i3,
     &       ' X =',1p,3e11.3)

3006  format(1x,'Plot',i3,' Stress : Elmt  =',i8,' No. =',i3)

3007  format(1x,'Plot',i3,' Contact: Slave =',i8,' DOF =',i3)

3008  format(1x,'Plot',i3,' React. : Node  =',i8,' DOF =',i3,
     &       ' X =',1p,3e11.3)

3009  format(1x,'Plot',i3,' Energy : Comp. =',i8,' Type=',i3)

3010  format(1x,'Plot',i3,' Arclen : Node  =',i8,' DOF =',i3)

3011  format(1x,'Plot',i3,' Sums   : DOF   =',i8,' DIR =',i3,
     &      ' X = ',1p,1e12.5,' TOL = ',1p,1e10.3)

3012  format(1x,'Plot',i3,' Matl.  : Model =',i8,' Type=',i3)

3013  format(1x,'Plot',i3,' Users  : Elmt  =',i8,' No. =',i3)

3014  format(1x,'Plot',i3,' R sum  : Dof   =',i8,' Node=',i8,
     &        ' Node=',i8)

      end
