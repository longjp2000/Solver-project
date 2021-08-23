c$Id: vblke.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine vblke(nr,ns,nt,x,ix,ni,ne,nf,ndm,nen1,mat,ntyp,
     &                 dlayer,ilr)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Generate a block of 3-d 8-node brick elements

c      Inputs:
c         nr        - Number elements in 1-local coordinate dir.
c         ns        - Number elements in 2-local coordinate dir.
c         nt        - Number elements in 3-local coordinate dir.
c         ni        - Initial node number for block
c         ne        - Initial element number for block
c         nf        - Final   element number for block
c         ndm       - Spatial dimension of coordinate array
c         nen1      - Dimension of ix array
c         mat       - Material set number for block
c         ntyp      - Element type for generations
c                     10:  8-node hexahedron  elements
c                     11:  4-node tetrahedron elements
c                     12: 27-node hexahedron  elements
c                     13: 10-node tetrahedron elements
c                     14: 20-node hexahedron  elements
c                     15: 11-node tetrahedron elements

c      Outputs:
c         x(ndm,*)  - Nodal coordinates 
c         ix(*)     - Element nodal connection list for block
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat2.h'
      include  'iofile.h'
      include  'trdata.h'

      integer   ni,nf,ne,nn,ndm,nen1,mat,ma,ntyp
      integer   nr,ns,nt,nrs,i,j,k,l,m,n,dlayer, ninc, nodi
      integer   ix(nen1,*),iq(27),it(10),ilr(*),nd(27),itl(4,6)
      integer   itq(10,6)
      real*8    x(ndm,*)

      save

      data      itl / 1,2,4,5, 2,3,4,8, 2,4,5,8,
     &                2,6,3,8, 3,6,7,8, 5,6,2,8 /

      data      itq / 1,2,4,5,  9, 21, 12, 17, 25, 23,
     &                2,3,4,8, 10, 11, 21, 27, 26, 20,
     &                2,4,5,8, 21, 23, 25, 27, 20, 16,
     &                2,6,3,8, 18, 24, 10, 27, 22, 26,
     &                3,6,7,8, 24, 14, 19, 26, 22, 15,
     &                5,6,2,8, 13, 18, 25, 16, 22, 27 /

c     Check generation order

      do i = 1,10
        it(i) = i
      end do ! i
      do i = 1,27
        iq(i) = i
      end do ! i
      if    (ntyp.eq.12) then
        nn = 27
      elseif(ntyp.eq.14) then
        nn = 20
      else
        nn = 8
      endif
      if(trdet.lt.0.0d0) then
        do i = 1,4
          iq(i+4) = i
          iq(i  ) = i+4
        end do ! i
        if(ntyp.eq.12 .or. ntyp.eq.14) then
          do i = 9,12
            iq(i+4) = i
            iq(i  ) = i+4
          end do ! i
          iq(21) = 22
          iq(22) = 21
        endif
        i     = it(2)
        it(2) = it(3)
        it(3) = i
        if(ntyp.eq.13 .or. ntyp.eq.15) then
          i      = it( 5)
          it( 5) = it( 7)
          it( 7) = i
          i      = it( 9)
          it( 9) = it(10)
          it(10) = i
        endif
      endif

      nrs = nr*ns
      if(ntyp.le.11) then
        ninc  =  1
        nd(1) = -1
        nd(2) =  0
        nd(3) =      nr
        nd(4) = -1 + nr
        nd(5) = -1      + nrs
        nd(6) =           nrs
        nd(7) =      nr + nrs
        nd(8) = -1 + nr + nrs
      elseif(ntyp.le.14) then
        ninc  =  2
        nd( 1) = -1
        nd( 2) =  1
        nd( 3) =  1 + nr*2
        nd( 4) = -1 + nr*2
        nd( 5) = -1        + nrs*2
        nd( 6) =  1        + nrs*2
        nd( 7) =  1 + nr*2 + nrs*2
        nd( 8) = -1 + nr*2 + nrs*2
        nd( 9) =  0
        nd(10) =  1 + nr
        nd(11) =      nr*2
        nd(12) = -1 + nr
        nd(13) =             nrs*2
        nd(14) =  1 + nr   + nrs*2
        nd(15) =      nr*2 + nrs*2
        nd(16) = -1 + nr   + nrs*2
        nd(17) = -1        + nrs
        nd(18) =  1        + nrs
        nd(19) =  1 + nr*2 + nrs
        nd(20) = -1 + nr*2 + nrs
        nd(21) =      nr
        nd(22) =      nr   + nrs*2
        nd(23) = -1 + nr   + nrs
        nd(24) =  1 + nr   + nrs
        nd(25) =             nrs
        nd(26) =      nr*2 + nrs
        nd(27) =      nr   + nrs
      endif

c     Compute element connections

      nodi = ni + nr*ns*nt - 1     ! Start number for internal nodes
      if(dlayer.ge.0) then
        ma = mat
      endif
      nf = ne - 1
      do k = 1,nt-1,ninc
        if(dlayer.eq.3) then
          ma = ilr(k)
        endif
        do j = 1,ns-1,ninc
          if(dlayer.eq.2) then
            ma = ilr(j)
          endif
          n = nr*(j-1 + ns*(k-1)) + ni
          do i = 1,nr-1,ninc
            if(dlayer.eq.1) then
              ma = ilr(i)
            endif
            n = n + 1

c           8-node hexahedral elements

            if(ntyp.eq.10) then
              nf = nf + 1
              ix(nen1,nf) = ma
              do m = 1,8
                ix(iq(m),nf) = n + nd(m)
              end do ! m

c           4-node tetrahedral elements

            elseif(ntyp.eq.11) then
              do l = 1,6
                nf = nf + 1
                ix(nen1,nf) = ma
                do m = 1,4
                  ix(it(m),nf) = n + nd(itl(m,l))
                end do ! m
              end do ! l

c           20 and 27-node hexahedral elements

            elseif(ntyp.eq.12 .or. ntyp.eq.14) then
              nf = nf + 1
              ix(nen1,nf) = ma
              do m = 1,nn
                ix(iq(m),nf) = n + nd(m)
              end do ! m

c           10-node tetrahedral elements

            elseif(ntyp.eq.13) then
              do l = 1,6
                nf = nf + 1
                ix(nen1,nf) = ma
                do m = 1,10
                  ix(it(m),nf) = n + nd(itq(m,l))
                end do ! m
              end do ! l

c           11-node tetrahedral elements

            elseif(ntyp.eq.15) then
              do l = 1,6
                nf = nf + 1
                ix(nen1,nf) = ma
                do m = 1,10
                  ix(it(m),nf) = n + nd(itq(m,l))
                end do ! m

c               Internal point

                nodi      = nodi + 1
                ix(11,nf) = nodi

                do m = 1,ndm
                  x(m,nodi) = 0.250d0*(x(m,ix(it( 5),nf))
     &                               + x(m,ix(it( 6),nf))
     &                               + x(m,ix(it( 7),nf))
     &                               + x(m,ix(it( 8),nf))
     &                               + x(m,ix(it( 9),nf))
     &                               + x(m,ix(it(10),nf)))
     &                      - 0.125d0*(x(m,ix(it( 1),nf))
     &                               + x(m,ix(it( 2),nf))
     &                               + x(m,ix(it( 3),nf))
     &                               + x(m,ix(it( 4),nf)))
                end do ! m

              end do ! l

            endif

            n = n + ninc - 1

          end do ! i
        end do ! j
      end do ! k

      end
