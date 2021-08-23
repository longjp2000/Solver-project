c$Id: sblke.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine sblke(nr,ns,x,ix,ni,ne,n,ndm,nen1,nodinc,ntyp,nm,mat,
     &                 dlayer,ilr,ctype)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use constants from 'pconstant.h'                   14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: elements for 2-d problems

c         nm < 4 : Generates a line of elements
c           ns = 1     2-node line
c           ns = 2     3-node line
c           ns = <10  ns-node line
c         nm > 3 : Generates a block of elements
c           ntyp = 0   4-node quadrilaterals
c           ntyp = 1   3-node triangles - diags ll to ur
c           ntyp = 2   3-node triangles - diags ul to lr
c           ntyp = 3   3-node triangles - diags
c           ntyp = 4   3-node triangles - diags
c           ntyp = 5   3-node triangles - diags union jack - 1
c           ntyp = 6   3-node triangles - diags union jack - 2
c           ntyp = 7   6-node triangles - diags ll to ur
c           ntyp = 8   8-node quadrilaterals
c           ntyp = 9   9-node quadrilaterals
c           ntyp = 16 16-node quadrilaterals

c           ntyp =-1   3-node triangles - crossed pattern
c           ntyp =-7   7-node triangles - diags ll to ur (bubble node)

c      Inputs:
c         nr        - Number nodes in 1-local coordinate dir.
c         ns        - Number nodes in 2-local coordinate dir.
c         ni        - Initial node number for block
c         ne        - Initial element number for block
c         n         - Number of previous last node on block
c         ndm       - Spatial dimension of mesh
c         nen1      - Dimension of ix array
c         nodinc    - Increment to node numbers at end of each line
c         ntyp      - Block type
c         nm        - Number of nodes on block
c         mat       - Material set number for block
c         ctype     - Type of block coordinates

c      Outputs:
c         n         - Number of last node on block (added by 7-node tris).
c         x(ndm,*)  - Nodal coordinates for block
c         ix(*)     - Element nodal connection list for block
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'cdat2.h'
      include   'iofile.h'
      include   'pconstant.h'
      include   'trdata.h'

      logical    ityp, pcomp, bernfl
      character  ctype*15
      integer    nr,ns,ni,ne,ndm,nen1,nodinc,ntyp,mat,ma,dlayer
      integer    i,j,n,me,np,nm,nn,n8,inc,ix(nen1,*),iq(9),it(6),ilr(*)
      integer    ic(16),icr(16)
      real*8     rr,sn,cn, x(ndm,*)

      save

      data       icr/ 4, 3, 2, 1,10, 9, 8, 7, 6, 5, 12,11,16,15,14,13/

c     Set Bernstein flag

      bernfl = pcomp(ctype,'bern',4)

c     Generate elements

      if(dlayer.ge.0) then
        ma = mat
      endif
      if(ne.gt.0) then
        do i = 1,9
          iq(i) = i
        end do ! i
        do i = 1,16
          ic(i) = i
        end do ! i
        do i = 1,6
          it(i) = i
        end do ! i
        if(trdet.lt.0.0d0) then
          if(ntyp.eq.16) then
            do i = 1,16
              ic(i) = icr(i)
            end do ! i
          else
            iq(1) = 4
            iq(2) = 3
            iq(3) = 2
            iq(4) = 1
            iq(5) = 7
            iq(7) = 5
            it(1) = 2
            it(2) = 1
            it(5) = 6
            it(6) = 5
          endif
        endif
        me = ne - 1

c       Line generations

        if(bernfl .or. nm.lt.4) then

c         Check global mapping for transformation

          if(trdet.lt.0.0d0) then
            iq(1) = 2
            iq(2) = 1
          else
            iq(1) = 1
            iq(2) = 2
          endif

c         Generate elements

          inc = max(1,min(9,ns))
          nn  = ni
          do i = 1,nr-1,inc
            if(dlayer.eq.1) then
              ma = ilr(inc*i-inc+1)
            endif
            nn = nn + 1
            me = me + 1
            ix(nen1,me)  = ma
            ix(iq(1),me) = nn - 1
            if(inc.gt.1) then
              do j = 3,inc+1
                ix(j,me) = nn
                nn       = nn + 1
              end do ! j
            endif
            ix(iq(2),me) = nn
          end do ! i

c       Block generations

        elseif(ntyp.ge.0 .or. ntyp .eq. -7) then
          inc = 1
          if(abs(ntyp).ge.7) inc = 2
          if(ntyp.eq.16) inc = 3
          do j = 1,ns-1,inc
            if(dlayer.eq.2) then
              ma = ilr(inc*j-inc+1)
            endif
            if(ntyp.eq.8) then
              nn = (nr + nodinc*2 + (nr+1)/2)*(j-1)/2 + ni
              n8 =  nn
            else
              nn = (nr + nodinc)*(j-1) + ni
            endif
            do i = 1,nr-1,inc
              if(dlayer.eq.1) then
                ma = ilr(inc*i-inc+1)
              endif
              nn = nn + 1
              me = me + 1
              ix(nen1,me) = ma
              if(ntyp.eq.0) then
                if(ndm.eq.1) then
                  ix(1,me)     = nn - 1
                  ix(2,me)     = nn
                else
                  ix(iq(1),me) = nn - 1
                  ix(iq(2),me) = nn
                  ix(iq(3),me) = nn + nr + nodinc
                  ix(iq(4),me) = nn + nr - 1 + nodinc
                endif
              elseif(abs(ntyp).eq.7) then
                ix(it(1),me) = nn - 1
                ix(it(2),me) = nn + 1
                ix(it(3),me) = 2*(nr+nodinc) + nn + 1
                ix(it(4),me) = nn
                ix(it(5),me) = nr+nodinc + nn + 1
                ix(it(6),me) = nr+nodinc + nn
                if(ntyp.eq.-7) then
                  ix(7,me) = n
                  x(1,n)   = (x(1,ix(it(4),me))
     &                     +  x(1,ix(it(5),me))
     &                     +  x(1,ix(it(6),me)))*four9
     &                     - (x(1,ix(it(1),me))
     &                     +  x(1,ix(it(2),me))
     &                     +  x(1,ix(it(3),me)))*one9
                  x(2,n)   = (x(2,ix(it(4),me))
     &                     +  x(2,ix(it(5),me))
     &                     +  x(2,ix(it(6),me)))*four9
     &                     - (x(2,ix(it(1),me))
     &                     +  x(2,ix(it(2),me))
     &                     +  x(2,ix(it(3),me)))*one9
                  if(ndm.ge.2 .and. pcomp(ctype,'pola',4)) then
                    call pdegree(x(2,n), sn,cn)
                    rr     = x(1,n)
                    x(1,n) = x0(1) + rr*cn
                    x(2,n) = x0(2) + rr*sn
                  endif
                  n        = n + 1
                endif
                me           = me + 1
                ix(it(1),me) = nn - 1
                ix(it(2),me) = 2*(nr+nodinc) + nn + 1
                ix(it(3),me) = 2*(nr+nodinc) + nn - 1
                ix(it(4),me) = nr+nodinc + nn
                ix(it(5),me) = 2*(nr+nodinc) + nn
                ix(it(6),me) = nr+nodinc + nn - 1
                ix(nen1,me)  = ma
                if(ntyp.eq.-7) then
                  ix(7,me) = n
                  x(1,n)   = (x(1,ix(it(4),me))
     &                     +  x(1,ix(it(5),me))
     &                     +  x(1,ix(it(6),me)))*four9
     &                     - (x(1,ix(it(1),me))
     &                     +  x(1,ix(it(2),me))
     &                     +  x(1,ix(it(3),me)))*one9
                  x(2,n)   = (x(2,ix(it(4),me))
     &                     +  x(2,ix(it(5),me))
     &                     +  x(2,ix(it(6),me)))*four9
     &                     - (x(2,ix(it(1),me))
     &                     +  x(2,ix(it(2),me))
     &                     +  x(2,ix(it(3),me)))*one9
                  if(ndm.ge.2 .and. pcomp(ctype,'pola',4)) then
                    call pdegree(x(2,n), sn,cn)
                    rr     = x(1,n)
                    x(1,n) = x0(1) + rr*cn
                    x(2,n) = x0(2) + rr*sn
                  endif
                  n        = n + 1
                endif
                nn = nn + 1
              elseif(ntyp.eq.8) then
                ix(iq(1),me) = nn - 1
                ix(iq(2),me) = nn + 1
                ix(iq(3),me) = nr + nodinc + (nr+1)/2 + nodinc + nn + 1
                ix(iq(4),me) = nr + nodinc + (nr+1)/2 + nodinc + nn - 1
                ix(iq(5),me) = nn
                ix(iq(6),me) = nr + nodinc +  n8 + 1
                ix(iq(7),me) = nr + nodinc + (nr+1)/2 + nodinc + nn
                ix(iq(8),me) = nr + nodinc +  n8
                nn = nn + 1
                n8 = n8 + 1
              elseif(ntyp.eq.9) then
                ix(iq(1),me) = nn - 1
                ix(iq(2),me) = nn + 1
                ix(iq(3),me) = 2*(nr+nodinc) + nn + 1
                ix(iq(4),me) = 2*(nr+nodinc) + nn - 1
                ix(iq(5),me) = nn
                ix(iq(6),me) = nr+nodinc + nn + 1
                ix(iq(7),me) = 2*(nr+nodinc) + nn
                ix(iq(8),me) = nr+nodinc + nn - 1
                ix(iq(9),me) = nr+nodinc + nn
                nn = nn + 1
              elseif(ntyp.eq.16) then
                ix(ic( 1),me) = nn - 1
                ix(ic( 2),me) = nn + 2
                ix(ic( 3),me) = 3*(nr+nodinc) + nn + 2
                ix(ic( 4),me) = 3*(nr+nodinc) + nn - 1
                ix(ic( 5),me) = nn
                ix(ic( 6),me) = nn + 1
                ix(ic( 7),me) =    nr+nodinc  + nn + 2
                ix(ic( 8),me) = 2*(nr+nodinc) + nn + 2
                ix(ic( 9),me) = 3*(nr+nodinc) + nn + 1
                ix(ic(10),me) = 3*(nr+nodinc) + nn
                ix(ic(11),me) = 2*(nr+nodinc) + nn - 1
                ix(ic(12),me) =    nr+nodinc  + nn - 1
                ix(ic(13),me) =    nr+nodinc  + nn
                ix(ic(14),me) =    nr+nodinc  + nn + 1
                ix(ic(15),me) = 2*(nr+nodinc) + nn + 1
                ix(ic(16),me) = 2*(nr+nodinc) + nn
                nn = nn + 2
              else
                ityp = (ntyp.eq.1)  .or.
     &                 (ntyp.eq.3.and.mod(j,2)  .eq.1) .or.
     &                 (ntyp.eq.4.and.mod(j,2)  .eq.0) .or.
     &                 (ntyp.eq.5.and.mod(i+j,2).eq.0) .or.
     &                 (ntyp.eq.6.and.mod(i+j,2).eq.1)
                if(ityp) then
                  ix(it(1),me) = nn - 1
                  ix(it(2),me) = nn + nr + nodinc
                  ix(it(3),me) = nn + nr + nodinc - 1
                  me = me + 1
                  ix(it(1),me) = nn + nr + nodinc
                  ix(it(2),me) = nn - 1
                  ix(it(3),me) = nn
                  ix(nen1,me)  = ma
                else
                  ix(it(1),me) = nn
                  ix(it(2),me) = nn + nr + nodinc - 1
                  ix(it(3),me) = nn - 1
                  me = me + 1
                  ix(it(1),me) = nn + nr + nodinc - 1
                  ix(it(2),me) = nn
                  ix(it(3),me) = nn + nr + nodinc
                  ix(nen1,me)  = ma
                endif
              endif
            end do ! i
          end do ! j
        elseif(ntyp.eq. -1) then
          do j = 1,ns-1
            if(dlayer.eq.2) then
              ma = ilr(j)
            endif
            n  = (2*nr+nodinc-1)*(j-1) + ni
            nn = n + 2*nr + nodinc - 1
            inc = 2
            if(j.eq.ns-1) inc = 1
            do i = 1,nr-1
              if(dlayer.eq.1) then
                ma = ilr(i)
              endif
              n            = n  + 1
              np           = nn + inc
              me           = me + 1
              ix(nen1,me)  = ma
              ix(it(1),me) = n - 1
              ix(it(2),me) = n + 1
              ix(it(3),me) = n
              me           = me + 1
              ix(nen1,me)  = ma
              ix(it(1),me) = n + 1
              ix(it(2),me) = np
              ix(it(3),me) = n
              me           = me + 1
              ix(nen1,me)  = ma
              ix(it(1),me) = np
              ix(it(2),me) = nn
              ix(it(3),me) = n
              me           = me + 1
              ix(nen1,me)  = ma
              ix(it(1),me) = nn
              ix(it(2),me) = n - 1
              ix(it(3),me) = n
              n            = n  + 1
              nn           = nn + inc
            end do ! i
          end do ! j
        endif
      endif

c     Set final element number

      n = me

      end
