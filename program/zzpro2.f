c$Id: zzpro2.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine zzpro2(ix,ib,ip,iq, x, st, sp,
     &                  ndm,nen,nen1,numnp,numel,numst)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Zienkiewicz-Zhu Projections

c      Inputs:

c      Outputs:
c-----[--+---------+---------+---------+---------+---------+---------+-]


      implicit   none

      include   'eldata.h'
      include   'strnum.h'
      include   'pconstant.h'
      include   'pointer.h'
      include   'comblk.h'
      include   'zzcom1.h'

      logical    fa
      integer    ndm, nen,nen1,numnp,numel,numst
      integer    i,ii, j,jj, k,kk, l, nn, nps

      real*8     cd,press

      integer    ix(nen1,numel), ib(numnp), ip(numnp),iq(*)
      real*8     x(ndm,*),st(numnp,numst), sp(numnp,*)
      real*8     hh(10,10),aa(30,10),bb(30,10),poly(10),t(8), sig(9)
      real*8     tmp(10)

      save

      data       fa /.false./

c     Set patch size parameters

      do i = 1,numst
        do nn = 1,numnp
          st(nn,i) = 0.0d0
        end do ! nn
      end do ! i

c     5.) Build patches

      do nn = 1,ip(numnp)
        iq(nn) = 0
      end do ! nn

      do nn = 1,numel
        do i = 1,nen
          ii = ix(i,nn)
          if(ib(ii).eq.0) then
            if(ii.gt.1) then
              jj = ip(ii-1) + 1
            else
              jj = 1
            end if
 300        if(iq(jj).eq.0 .or. iq(jj).eq.nn) then
              iq(jj) = nn
            else
              jj = jj + 1
              go to 300
            end if
          end if
        end do ! i
      end do ! nn

c     6.) Loop over patches to compute projected quantities

      jj  = 0
      do nn = 1,numnp
        if(ib(nn).eq.0) then

          call pzero(hh,100)
          call pzero(aa,300)
          call pzero(bb,300)
          do i = 1,ndm
            xnodz(i) = x(i,nn)
          end do ! i

c         Form stresses for patch

          do j = jj+1,ip(nn)

c           Get contribution from each element

            if(iq(j).gt.0) then

c             Set default projection quantities and initialize

              nums = numst
              do i = 1,10
                do k = 1,10
                  ek(k,i) = 0.0d0
                end do ! k
                do k = 1,nums
                  est(k,i) = 0.0d0
                end do ! k
              end do ! i

c             Compute projections

              call formfe(np(40),np(26),np(26),np(26),fa,fa,fa,fa,25,
     &                    iq(j),iq(j),1)

              nps = 3

              if(ndm.eq.2) then
                if(nel.eq.6 .or. nel.eq.8 .or. nel.eq.9) then
                  nps = 6
                elseif(nel.eq.16) then
                  nps = 10
                endif
              endif

c             Assemble least square patch arrays

              do i = 1,nps
                do k = 1,nps
                  hh(k,i) = hh(k,i) + ek(k,i)
                end do ! k
                do k = 1,nums
                  aa(k,i) = aa(k,i) + est(k,i)
                end do ! k
              end do ! i

            endif

          end do ! j

c         Solve

          call zinvert(hh,nps,10,tmp)

          do i = 1,nums
            do j = 1,nps
              t(j) = 0.0d0
            end do ! j
            do k = 1,nps
              do j = 1,nps
                t(j) = t(j) + hh(j,k)*aa(i,k)
              end do ! j
            end do ! k
            do j = 1,nps
              bb(i,j) = bb(i,j) + t(j)
            end do ! j
          end do ! i

c         Assemble nodal stresses

          do j = jj+1,ip(nn)
            if(iq(j).gt.0) then
              do i = 1,nen
                kk = ix(i,iq(j))
                if(kk.gt.0) then
                  if(ib(kk).gt.0) then

c                   Increment ib(kk) to count number of nodal projections

                    ib(kk) = ib(kk) + 1

c                   Compute polynomial

                    poly(1) = 1.d0
                    do k = 1,ndm
                      poly(k+1) = x(k,kk) - x(k,nn)
                    end do ! k

                    if(ndm.eq.2) then
                      if(nps.ge.6) then
                        poly( 4) = poly(2)*poly(2)
                        poly( 5) = poly(2)*poly(3)
                        poly( 6) = poly(3)*poly(3)
                      endif
                      if(nps.eq.10) then
                        poly( 7) = poly(4)*poly(2)
                        poly( 8) = poly(5)*poly(2)
                        poly( 9) = poly(5)*poly(3)
                        poly(10) = poly(6)*poly(3)
                      endif
                    endif

c                   Project boundary nodes

                    do k = 1,nums
                      do l = 1,nps
                        st(kk,k) = st(kk,k) + bb(k,l)*poly(l)
                      end do ! l
                    end do ! k
                  end if
                end if
              end do ! i
            end if
          end do ! j

c         Assemble interior node patch value

          do i = 1,nums
            st(nn,i) = bb(i,1)
          end do ! i

        end if

        jj = ip(nn)

      end do ! nn

c     7.) Boundary averages

      do nn = 1,numnp
        if(ib(nn).gt.1) then

          cd = 1.d0/dble(ib(nn)-1)
          do i = 1,nums
            st(nn,i) = st(nn,i)*cd
          end do ! i
        end if

c       Three-dimensional Principal Stresses

        if(ndm.eq.3 .or. istp.eq.8) then
          sig(1)   = st(nn,1)
          sig(2)   = st(nn,2)
          sig(3)   = st(nn,3)
          sig(4)   = st(nn,4)
          sig(5)   = st(nn,5)
          sig(6)   = st(nn,6)
          call pstr3d(sig,sig(7))
          sp(nn,1) = sig(7)
          sp(nn,2) = sig(8)
          sp(nn,3) = sig(9)

c       Two-dimensional Principal Stresses

        elseif(ndm.eq.2) then
          sig(1)   = st(nn,1)
          sig(2)   = st(nn,2)
          sig(3)   = st(nn,3)
          sig(4)   = st(nn,4)
          call pstr2d(sig,sig(7))
          sp(nn,1) = sig(7)
          sp(nn,2) = sig(8)
          sp(nn,3) = (sig(7)-sig(8))*0.5d0
          sp(nn,4) = sig(9)
        endif

c       Compute mean stress and mises stress

        press    = (sp(nn,1) + sp(nn,2) + sp(nn,3))*one3
        sp(nn,5) = press
        sp(nn,6) = sqrt(1.5d0*((sp(nn,1) - press)**2
     &                       + (sp(nn,2) - press)**2
     &                       + (sp(nn,3) - press)**2))
        sp(nn,7) =        one3*(sp(nn,1) - press)**3
     &                        *(sp(nn,2) - press)**3
     &                        *(sp(nn,3) - press)**3
        if(sp(nn,7).lt.0.0d0) then
          sp(nn,7) = -(abs(sp(nn,7))**one3)
        else
          sp(nn,7) =  (abs(sp(nn,7))**one3)
        endif

      end do ! nn

      end
