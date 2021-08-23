c$Id: pelink.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pelink(id,x,ndm,ndf,numnp,neq,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Perform link of degree of freedom based on input
c               of edge coordinate values and direction.

c      Inputs:
c         x(ndm,*)  - Nodal coordinates of mesh
c         ndm       - Spatial dimension of mesh
c         ndf       - Number dof/node
c         numnp     - Number of nodes in mesh
c         prt       - Output results if true

c      Outputs:
c         id(ndf,*) - Modified equation numbers from links
c         neq       - Number of equations after links
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdat2.h'
      include  'comfil.h'
      include  'conval.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'ioincl.h'
      include  'part0.h'
      include  'trdata.h'

      logical   prt,lsave,errck,tinput,vinput,pcomp, oprt,prth, equalfl
      character text*15, fnamr*132,fext*4, type*4
      integer   ndm,ndf,numnp,neq, iosfil, i,ii,i1,i2, j, m1,m2,nmax
      integer   id(ndf,*),idl(12)
      real*8    gap, x1, x2, td(16),x(ndm,*)

      save

c     Routine to link degrees of freedom together

      gap   =  1.d-04
      fnamr =  fsav
      fext  =  'eln'
      call addext(fnamr,fext,128,4)
      call opnfil(fext,fnamr,-2,ios,lsave)
      if(lsave) then
        iosfil = ior
        ior    = ios

        oprt = prt
        do i = 0,36
          do j = 1,26
            vvsave(j,i) = vvv(j,i)
          end do ! j
        end do ! i
        do i = 1,3
          x0sav(i) = x0(i)
        end do ! i
        read(ior,1000) type,fincld(isf),irecrd(isf),prt,prth
        read(ior,1001) vvv
        read(ior,1001) tr,xr,trdet,x0

10      errck =  tinput(text,1,td(2),2+ndf)
        if(pcomp(text,'gap',3)) then
          gap   = td(2)
          if(prt) then
            write(iow,2001) gap
            if(iosfil.lt.0) then
              write(*,2001) gap
            endif
          endif
          go to 10
        else
          errck = vinput(text,15,td(1),1)
          i1  = td(1)
        endif
        i1  = min(i1,ndm)
        x1  = td(2)
        x2  = td(3)
        do i = 1,ndf
          idl(i) = td(i+3)
        end do ! i
        if(i1.eq.0) then
          close(ior)
          ior = iosfil

          prt = oprt
          do i = 0,36
            do j = 1,26
              vvv(j,i) = vvsave(j,i)
            end do ! j
          end do ! i
          do i = 1,3
            x0(i) = x0sav(i)
          end do ! i

          return
        endif

c       Look for matching coordinates in i1-direction

        do m1 = 1,numnp-1
          if(abs(x(i1,m1)-x1) .lt. gap ) then
            do m2 = m1+1,numnp
              if( abs(x(i1,m2)-x2) .lt. gap ) then

c               Check that all other values are same

                if(abs(x1-x2) .gt. gap) then
                  equalfl = .false.
                  do i2 = 1,ndm
                    if(i2.ne.i1) then
                      if(abs(x(i2,m1)-x(i2,m2)) .gt. gap ) go to 20
                    endif
                  end do ! i2
                else
                  equalfl = .true.
                endif

c               Match found print result

                if(prt) then
                  write(iow,2000) m1,m2,(idl(i),i=1,ndf)
                  if(iosfil.lt.0) then
                    write(*,2000) m1,m2,(idl(i),i=1,ndf)
                  endif
                endif
                do j = 1,ndf
                  if(ndfp(j).eq.npart .and. idl(j).eq.0) then
                    if(id(j,m1).gt.0 .and. id(j,m2).gt.0) then

c                     Select node to renumber dof

                      if(id(j,m1).lt.id(j,m2)) then
                        nmax     = id(j,m2)
                        id(j,m2) = id(j,m1)
                      else
                        nmax     = id(j,m1)
                        id(j,m1) = id(j,m2)
                      endif
                      do ii = 1,numnp
                        if(id(j,ii).eq.nmax) then
                          id(j,ii) = id(j,m1)
                        end if
                      end do ! ii

c                     Loop through all nodes to reduce equation numbers

                      errck = .false.
                      do i = 1,ndf
                        if(ndfp(i).eq.npart) then
                          do ii = 1,numnp
                            if(id(i,ii).gt.nmax) then
                              id(i,ii) = id(i,ii) - 1
                              errck    = .true.
                            endif
                          end do ! ii
                        endif
                      end do ! i
                      if(errck) neq = neq - 1

                    endif

                  end if
                end do ! j
              endif
20            continue
            end do ! m2
            if(equalfl) go to 10
          endif
        end do ! m1
        go to 10
      else
        write(iow,3000)
        write(ilg,3000)
        call plstop()
      endif

c     Formats

1000  format(a4,2x,a12,i8,2l5)
1001  format(1p,4e20.12)

2000  format(5x,'Link node',i8,' to',i8,' DOF = 0 to link:',12i4)
2001  format(/5x,'Search gap =',1p,1e12.4/)

3000  format(5x,'*ERROR* PELINK: Link file does not exist')

      end
