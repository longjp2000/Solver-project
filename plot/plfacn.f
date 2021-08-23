c$Id: plfacn.f,v 1.1 2006/11/20 20:33:12 rlt Exp $
      subroutine plfacn(ix,ia,nen,numel,nface,ie,nie)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Determines which exterior faces are directed toward
c               view point.

c      Inputs:
c         ix(nen1,*)- Element nodal connection lists
c         ia(*)     - Active element plots based on materials
c         nen       - Number nodes/element
c         numel     - Number of elements/faces
c         ie(nie,*) - Material set assembly data
c         nie       - Dimension of ie array

c      Outputs:
c         nface     - Number of faces
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'pbody.h'
      include  'pdata5.h'
      include  'pdata6.h'
      include  'plclip.h'
      include  'sdata.h'

      include  'pointer.h'
      include  'comblk.h'

      logical   lclip,addfac
      integer   nen,numel,nface,nie, i,j,m,n, iel,iiel, ien,nel,pstyp
      integer   ix(nen1,numel), ia(*), ie(nie,*), iq(4,6), it(3,4)

      save

c     8-node brick faces

      data iq/3,2,1,4, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 5,6,7,8/

c     4-node tet faces

      data it/1,2,4, 2,3,4, 3,1,4, 1,3,2/

c     Compute location of boundary faces

      nface = 0
      do n = 1,numel
        if(ix(nen1-1,n).ge.nreg1 .and. ix(nen1-1,n).le.nreg2 .and.
     &     ix(nen1-1,n).ge.0 .and. ia(n).ge.0) then
         pstyp = ie(1,ix(nen1,n))
         if(pstyp.gt.0) then
           iel = ie(nie-1,ix(nen1,n))
           do j = nen,1,-1
             if(ix(j,n).gt.0) then
               nel = j
               exit
             endif
           end do ! j

c          Get plot type

           call plftyp(pstyp,nel,iel)

           if(iel.gt.0) then
             iiel = inord(iel)
           else
             iiel = exord(-iel)
           endif

c          6-node triangle

           if(iiel.eq.7) then
             ien = 3
           else
             ien = nen
           endif

c          No face if inord < 0

           if     (iiel.lt.0) then

c          Set for tetrahedral element faces

           elseif (iiel .eq. 9 .or. iiel .eq. 15) then

            if( lclip(ix(1,n),4,hr(np(43)),ndm) ) then
              do m = 1,4
                addfac = .true.
                do j = 1,3
                  i = ix(it(j,m),n) - 1
                  if(mr(np(89)+i).eq.0) then
                    addfac = .false.
                  endif
                end do ! j
                if(addfac) then
                  nface = nface + 1
                endif
              end do ! m
            end if

c          Set for brick element faces

           elseif (iiel .gt. 10 ) then

            if( lclip(ix(1,n),8,hr(np(43)),ndm) ) then
              do m = 1,6
                addfac = .true.
                do j = 1,4
                  i = ix(iq(j,m),n) - 1
                  if(mr(np(89)+i).eq.0) then
                    addfac = .false.
                  endif
                end do ! j
                if(addfac) then
                  nface = nface + 1
                endif
              end do ! m
            end if

c          Set space for line elements

           elseif( iiel.gt.0 .and. iiel.le.3 ) then

             if( lclip(ix(1,n),2,hr(np(43)),ndm) ) then
               nface = nface + 1
             endif
  
c          Set space for top and bottom shell faces

           elseif( lclip(ix(1,n),min(4,ien),hr(np(43)),ndm) ) then

            nface = nface + 2

           end if

          end if
        end if ! pty
      end do ! n

      end
