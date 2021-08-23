c$Id: umacr5.f,v 1.1 2006/11/21 16:46:48 rlt Exp $
      subroutine umacr5(lct,ctl,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c     Authors: D.S. Bindel, S. Govindjee 6/21/2005
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Serial ARPACK interface to FEAP

c      Inputs:
c         lct        - Command character parameters
c                      'symm'etric: (K - lambda M)*V = 0
c         ctl(3)     - Command numerical parameters
c                      1 = Number eigenpairs to compute
c                      2 = Maximum number of iterations
c                      3 = Tolerance for eigenvalues
c         prt        - Flag, output if true

c      Outputs:
c         hr(np(76)) - Eigenvalues
c         hr(np(77)) - Eigenvectors
c      Purpose:  User interface for adding solution command language
c                instructions.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'evdata.h'
      include  'iofile.h'
      include  'rdata.h'
      include  'umac1.h'
      include  'p_point.h'
      include  'part0.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   pcomp,prt, setval,palloc
      character lct*15
      real*8    ctl(3),stol

c     ARPACK STUFF BEGIN

      integer   maxitr, mmode

c     ARPACK STUFF END

      save

c     Set command word

      if(pcomp(uct,'mac5',4)) then      ! Usual    form
        uct = 'arpa'                    ! Specify 'arpa'ck
      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

      else

        mf     = min(neq,max(1,nint(ctl(1))))
        mq     = min(neq,max(20,2*mf))
        maxitr = max(300,nint(ctl(2)))
c       vneq   = neq
        if(ctl(3).ne.0.0d0) then
          stol = ctl(3)
        else
          stol = max(1.d-12,tol)
        endif

c       Symmetric eigensolve

        if(pcomp(lct,'symm',4) .or. pcomp(lct,'lump',4)
     &                         .or. pcomp(lct,'    ',4)) then

          if(pcomp(lct,'lump',4)) then
            mmode = 1
            setval = palloc(157,'USER7',neq,2)    ! Sqrt M
            point  = np(157)
            call masssqr(hr(point),hr(np(12+npart)), neq)  
          else
            point = np(1)
            mmode = 3
          endif

          if(ior.lt.0) then
            write(  *,2001) 'Symmetric',mf,mq,shift,maxitr,neq
          else
            write(iow,2001) 'Symmetric',mf,mq,shift,maxitr,neq
          endif

          setval = palloc( 76,'EVAL ', 2*mq,         2)  ! d
          setval = palloc( 77,'EVEC ', mq*neq ,      2)  ! v
          setval = palloc(158,'USER8', neq,          2)  ! resid
          setval = palloc(159,'USER9', mq*(mq+8),    2)  ! workl
          setval = palloc(160,'USER0', 3*neq,        2)  ! workd

          call arfeaps(mmode,      hr(np( 76)), hr(np( 77)),
     &                 hr(point)  ,hr(np(158)), hr(np(159)),
     &                 hr(np(160)),neq, mf, mq, maxitr, shift,
     &                 stol, prt)

          setval = palloc(160,'USER7',  0,2)    ! Delete work2
          setval = palloc(159,'USER7',  0,2)    ! Delete work1
          setval = palloc(158,'USER7',  0,2)    ! Delete resid
          if(mmode.eq.1) then
            setval = palloc(157,'USER7',  0,2)    ! Delete mass sqrt
          endif

        else
           write(*,*) '--> ARPACK: Please Specify Option: symm or iter'
           return
        end if

c       Delete common scratch arrays

        setval = palloc(158,'USER8', 0,2)  ! resid
        setval = palloc(159,'USER9', 0,2)  ! workl
        setval = palloc(160,'USER0', 0,2)  ! workd

      end if

c     Format

2001  format(/10x,a,' Eigenproblem',/,
     &        12x,'No. vecs requested  ',i6,/
     &        12x,'No. vecs being used ',i6,/
     &        12x,'Shift:              ',1p,1e15.8,/
     &        12x,'Maximum # iterations',i6,/
     &        12x,'Number eqns. (neq)  ',i6)

      end

      subroutine masssqr(msqrt,mass, neq)  

c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Compute reciprocal of square root of mass

c     Inputs:
c       mass(*)  - Lumped mass
c       neq      - Number of equations

c     Outputs:
c       msqrt(*) - Reciprocal of square root of mass
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      integer    neq,n
      real*8     msqrt(*), mass(*)

      do n = 1,neq
        if(mass(n).gt.0.0d0) then
          msqrt(n) = 1.d0/sqrt(mass(n))
        else
          msqrt(n) = 1.d0
        endif
      end do

      end
