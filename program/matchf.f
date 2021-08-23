c$Id: matchf.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      logical function matchf(factyp,nec, mi,ix)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Set matching interfaces

c     Inputs:
c       factyp - Type of interface
c       nec    - Check unumber
c       mi     - Position in list to check
c       ix(*)  - List of element nodes

c     Outputs:
c       matchf - .true. if a match occurs
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'facset.h'
      include   'iofile.h'

      integer    factyp,nec,mi,mj,mk, ix(*), me(4,6), jchk(3)

      save

      data   me / 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 1,2,3,4, 5,6,7,8/

      if(factyp.eq.1) then
        jchk(1) = ix(mi)
        matchf  = ichk(1).eq.jchk(1)
      elseif(factyp.eq.2) then
        mj      = mod(mi,nec) + 1
        jchk(1) = min(ix(mi),ix(mj))
        jchk(2) = max(ix(mi),ix(mj))
        matchf  = ichk(1).eq.jchk(1) .and. ichk(2).eq.jchk(2)

      elseif(factyp.eq.3) then
        mj      = mod(mi,nec) + 1
        mk      = mod(mj,nec) + 1
        jchk(1) = min(ix(mi),ix(mj),ix(mk))
        jchk(2) = max(ix(mi),ix(mj),ix(mk))
        jchk(3) = ix(mi) + ix(mj) + ix(mk) - ichk(1) - ichk(2)
        matchf  = ichk(1).eq.jchk(1) .and. ichk(2).eq.jchk(2)
     &                               .and. ichk(3).eq.jchk(3)
      elseif(factyp.eq.4) then
        jchk(1) = min(ix(me(1,mi)),ix(me(2,mi)),
     &                ix(me(3,mi)),ix(me(4,mi)))
        do mj = 1,4
          if(ix(me(1,mi)).eq.ichk(1)) then
            mk      = mod(mj+1,4) + 1
            jchk(2) = ix(me(mk,mi))
            exit
          endif
        end do ! mj
        matchf  = ichk(1).eq.jchk(1) .and. ichk(2).eq.jchk(2)
      else
        matchf = .false.
        write(ilg,3000) factyp
      endif

c     Format

3000  format(' *ERROR* in MATCHF: No facetyp =',i5)
      end
