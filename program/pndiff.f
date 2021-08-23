c$Id: pndiff.f,v 1.2 2006/12/08 00:12:35 rlt Exp $
      subroutine pndiff(ix,ul,st,k1,du,flg)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute tangent by numerical differentiation

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat1.h'
      include  'ddata.h'
      include  'debugs.h'
      include  'eldata.h'
      include  'eltran.h'
      include  'iofile.h'
      include  'part0.h'
      include  'pointer.h'
      include  'sdata.h'
      include  'tdata.h'
      include  'tdatb.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   flg
      integer   i,j,k,m, k1, ix(nen)
      real*8    du,sqrtdig,usize,denom,temp1,temp2,temp3,temp4,temp5
      real*8    ul(ndf,nen,*),st(nst,nst,2)

      save

c     Tolerance for numerical diffentiation

      data      sqrtdig /1.d-8/

c     Numerical differentiation to compute element tangent

      call pzero(st,nst*nst)

c     Compute current residual

      if(flg) then
        call pform(ul,hr(np(44)),hr(np(39)),mr(np(34)),
     &             hr(np(35)),hr(np(36)),mr(np(32)),hr(np(25)),
     &             mr(np(31)),hr(np(43)),mr(np(33)),mr(np(181)),
     &             hr(np(30)),hr(np(38)),mr(np(20+npart)),hr(np(40)),
     &             hr(np(42)),hr(np(26)),hr(np(26)),hr(np(26)),
     &             ndd,nie,ndf,ndm,nen1,nst,
     &             .false.,.false.,.false.,.true.,3,k1,k1,1)

c       Save element tangent

        do m = 1,nst
          fp(1) = np(36) + (m-1)*nst - 1
          do i = 1,nst
            st(i,m,2) = hr(fp(1)+i)
          end do ! i
        end do ! m
      endif

      usize = 0.0d0
      do i = 1,nen
        do j = 1,ndf
          usize = max(usize,abs(ul(j,i,1)))
        end do ! j
      end do ! i
      if(debug) then
        write(iow,*) ' Computed max U:',usize,cc1,cc2,c3,c4,c5,du
      endif

      if(du.ne.0.0d0) then
        usize = du
      elseif(usize.ne.0.0d0) then
        usize = sqrtdig*usize
      else
        usize = sqrtdig
      endif
      denom = 1.d0/usize
      denom = 0.5d0/usize

c     For each dof compute perturbed residual

      do i = 1,nen
        do j = 1,ndf
          temp1     = ul(j,i,1)
          temp2     = ul(j,i,2)
          temp3     = ul(j,i,3)
          temp4     = ul(j,i,4)
          temp5     = ul(j,i,5)
          ul(j,i,1) = temp1 + cc1*usize
          ul(j,i,2) = temp2 + cc2*usize
          if(np(42).ne.0) then
            ul(j,i,1) = temp1 + c3*usize
            if(noi.eq.5) then
              ul(j,i,4) = temp4 + c4*usize
              ul(j,i,5) = temp5 + c5*usize
            endif
          endif

c         Compute perturbed residual

          call elmlib(hr(np(25)),ul,hr(np(44)),ix,hr(np(39)),
     &                hr(np(36)),hr(np(35)),ndf,ndm,nst,iel,6)

c         Compute difference for tangent

          do k=1,nst
            st(k,ndf*(i-1)+j,1) = -hr(np(35)+k-1)
          end do

c         Set displacements back to original values

          ul(j,i,1)  = temp1
          ul(j,i,2)  = temp2
          ul(j,i,3)  = temp3
          ul(j,i,4)  = temp4
          ul(j,i,5)  = temp5

c         Compute perturbed displacements

          ul(j,i,1) = temp1 - cc1*usize
          ul(j,i,2) = temp2 - cc2*usize
          if(np(42).ne.0) then
            ul(j,i,1) = temp1 - c3*usize
            if(noi.eq.5) then
              ul(j,i,4) = temp4 - c4*usize
              ul(j,i,5) = temp5 - c5*usize
            endif
          endif

c         Compute perturbed residual

          call elmlib(hr(np(25)),ul,hr(np(44)),ix,hr(np(39)),
     &                hr(np(36)),hr(np(35)),ndf,ndm,nst,iel,6)

c         Compute difference for tangent

          do k=1,nst
            st(k,ndf*(i-1)+j,1) = (st(k,ndf*(i-1)+j,1)
     &                          +  hr(np(35)+k-1))*denom
          end do

c         Set displacements back to original values

          ul(j,i,1)  = temp1
          ul(j,i,2)  = temp2
          ul(j,i,3)  = temp3
          ul(j,i,4)  = temp4
          ul(j,i,5)  = temp5
        end do ! j
      end do ! i

c     Output numerical tangent

      if(flg) then
        call mprint(st(1,1,2),nst,nst,nst,'S_tangent')
        call mprint(st(1,1,1),nst,nst,nst,'S_numerical')

c       Compute and output difference in tangent

        do n = 1,nst
          do i = 1,nst
            st(i,n,1) = st(i,n,2) - st(i,n,1)
          end do ! i
        end do ! n

        call mprint(st,nst,nst,nst,'Difference: S_t - S_num')


c     Replace tangent by numerical tangent

      else
        do m = 1,nst
          j     = np(36) + (m-1)*nst - 1
          do i = 1,nst
            hr(j+i) = st(i,m,1)
          end do ! i
        end do ! m
      endif

      end
