c$Id: stfqtq2.f,v 1.1 2006/11/20 20:32:46 rlt Exp $
      subroutine stfqtq2 (cp0,cm,ch2,ch1,ch3,tanm,resv,xm)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronym: STiFfness for QTQ contact with 2 node edge state

c      Purpose: Compute tangent and residual

c      Inputs :
c         cp0(*)  - Contact pair control data
c         cm(*)   - Contact materials data storage
c         ch1(*)  - Contact history variables (old)
c         ch2(*)  - Contact history variables (current)
c         ch3(*)  - Contact history variables
c         xm(*,*) - Master node coordinates (current)

c      Outputs:
c         ch3(*)  - Contact history variables (current)
c         tanm(*) - Contact tangent array
c         resv(*) - Contact residual array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include   'c_0.h'
      include   'c_comnd.h'
      include   'c_keyh.h'
      include   'c_pair.h'

      integer   i,j,k,l,istgn,istgt,istgnold
      real*8    cp0(nr0,n0c3:*),cm(*),ch1(*),ch2(*),ch3(*)
      real*8    tanm(18,18),resv(3,4),xm(3,4),epsn,gn,a1(3)
      real*8    a11,nvec(3),a1oben(3),ninew(3),niablnew(3),epst
      real*8    ni(2),niabl(2),xi,new(3,3)
      real*8    delxi(3,3),delgn(3,3),ddgn(3,3,3,3),xabl(3)
      real*8    ttvec(3),ttnorm,mue,ftrial,diffxi(2)
      real*8    delttvec(3,3,3),delttrial(3,3,3),delta(3,3)

      save

      data      delta /1.d0, 3*0.0d0, 1.d0, 3*0.0d0, 1.d0 /

      call cdebug0 ('     stfqtq2',-1)

c     Get geometrical variables

      istgn    =  ch2(p1( 4))
      istgnold =  ch1(p1( 4))
      gn       =  ch3(p3(9))
      nvec(1)  =  ch3(p3(18))
      nvec(2)  =  ch3(p3(18)+1)
      nvec(3)  =  ch3(p3(18)+2)
      a1(1)    =  ch2(p1(19)  )
      a1(2)    =  ch2(p1(19)+1)
      a1(3)    =  ch2(p1(19)+2)
      a11      =  a1(1)*a1(1) + a1(2)*a1(2) + a1(3)*a1(3)
      xi       =  ch2(p1(24)  )
      ni(1)    =  1.d0 - xi
      ni(2)    =  xi
      niabl(1) = -1.d0
      niabl(2) =  1.d0

      do i=1,3
        xabl(i) = xm(i,1)*niabl(1) + xm(i,2)*niabl(2)
      enddo

c     Compute var. of xi

      ninew(1)    =  1.d0
      ninew(2)    = -ni(1)
      ninew(3)    = -ni(2)
      niablnew(1) =  0.d0
      niablnew(2) =  niabl(1)
      niablnew(3) =  niabl(2)

      do i = 1,3
        do j = 1,3
          delgn(i,j) =  nvec(i)*ninew(j)
          delxi(i,j) = (a1(i)*ninew(j) + gn*nvec(i)*niablnew(j))/a11
         enddo
      enddo

      do k = 1,3
        do j = 1,3
          new(j,k) = niablnew(j)*nvec(k)
        enddo
      enddo

      do i=1,3
        a1oben(i)=a1(i)/a11
      enddo

      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              ddgn(i,k,j,l)=-ninew(l)*(a1oben(j)*new(k,i))
     &                     - delxi(i,k)*niablnew(l)*nvec(j)
               enddo
            enddo
         enddo
      enddo

c     Get penalty for normal stiffness

      epsn = cp0(3,2)*ch3(p3(16))
      do i = 1,3
        do k = 1,3
          do j = 1,3
            do l = 1,3
              tanm(3*j-3+i,3*l-3+k) = epsn*(delgn(i,j)*delgn(k,l)
     &                              + gn*ddgn(i,j,k,l))
            end do
          end do
        end do
      end do
      do i = 1,3
        do j=1,3
          resv(i,j) = -epsn*gn*delgn(i,j)
        end do
      end do

      ch3(p3(14)) = cp0(3,2)*gn

      if (iffric.eq.1) then

        mue   = cm(1)
        istgt = ch2(p1 (3))
        epst  = cp0(4,2)*ch3(p3(16))

        if (istgt.eq.1) then

          if(istgnold.eq.1) then
            epst = cp0(4,2)*ch3(p3(16))
            write(*,*)'strange corner case'
            do i=1,3
              do j=1,3
                resv(i,j) = resv(i,j) - ch2(p1(26)+i-1)*epst*ninew(j)
              enddo
            enddo
          elseif((istgnold.eq.2).or.(istgnold.eq.3)) then

            if((ch2(p1( 1)).eq.ch1(p1( 1))) .and.
     &         (ch2(p1(25)).eq.ch1(p1(25)))) then

c             At same edge

              diffxi(1) = ch2(p1(10))
              epst      = cp0(4,2)*ch2(p1(26))

              do i = 1,3
                ttvec(i) = -epst*diffxi(1)*a1(i)
              enddo

              do i=1,3
                do j=1,3
                  do k=1,3
                    delttrial(i,j,k)=-(delta(i,k)*diffxi(1)*ninew(j)
     &                               + delxi(i,j)*a1(k))*epst
                  enddo
                enddo
              enddo

              ttnorm = sqrt(ttvec(1)**2 + ttvec(2)**2 + ttvec(3)**2)
              ftrial = ttnorm - mue*abs(epsn*gn)

              if (ftrial.gt.0) then

c               sliding

                do i=1,3
                  do j=1,3
                    do k=1,3
                      delttvec(i,j,k)= (epsn*delgn(i,j)*ttvec(k)
     &                               + abs(epsn*gn)*delttrial(i,j,k))
     &                               * mue/ttnorm
                      do l=1,3
                        delttvec(i,j,k) = delttvec(i,j,k)
     &                                  + ttvec(l)*delttrial(i,j,l)
     &                                  * abs(epsn*gn)*ttvec(k)
     &                                  * mue/(ttnorm*ttnorm*ttnorm)
                      enddo
                    enddo
                  enddo
                enddo
                do i=1,3
                  ttvec(i)=mue*ttvec(i)*abs(epsn*gn)/ttnorm
                enddo
              endif

              do i = 1,3
                do k=1,3
                  do j=1,3
                    do l=1,3
                      tanm(3*j-3+i,3*l-3+k) = tanm(3*j-3+i,3*l-3+k)
     &                                      - delttvec(i,j,k)*ninew(l)
                    end do
                  end do
                end do
              end do

              do i=1,3
                do j=1,3
                  resv(i,j)=resv(i,j) + ttvec(i)*ninew(j)
                enddo
              enddo

              ch3(p3( 5)) = ttvec(3)*ninew(1)
              ch3(p3(15)) = sqrt(ttvec(1)**2+ttvec(2)**2+ttvec(3)**2)

            else

              do i=1,3
                do j=1,3
                  resv(i,j)=resv(i,j) - ch2(p1(26)+i-1)*epst*ninew(j)
                enddo
              enddo

              ch3(p3( 5)) =-ch2(p1(26)+2)*epst*ninew(1)
              ch3(p3(15)) = sqrt(ch2(p1(26)  )**2
     &                         + ch2(p1(26)+1)**2
     &                         + ch2(p1(26)+2)**2)*epst

            endif
          elseif(istgnold.eq.4) then
            do i=1,3
              do j=1,3
                resv(i,j)=resv(i,j) - ch2(p1(26)+i-1)*epst*ninew(j)
              enddo
            enddo
            ch3(p3( 5)) =-ch2(p1(26)+2)*epst*ninew(1)
            ch3(p3(15)) = sqrt(ch2(p1(26)  )**2
     &                       + ch2(p1(26)+1)**2
     &                       + ch2(p1(26)+2)**2)*epst

          elseif(istgnold.eq.0) then

            ch3(p3( 5)) = 0.0d0
            ch3(p3(15)) = 0.0d0

          endif
        endif                  !istgt.eq.1
      endif                    !mue.ge.0.d0

      end
