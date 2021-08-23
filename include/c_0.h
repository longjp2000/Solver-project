
                               ! CONTACT PARAMETERS

      integer    c_ncc,c_ncs,c_ncel,c_lp1,c_lp3,c_lmv

      parameter (c_ncc = 12)   ! # of available contact commands
      parameter (c_ncs = 200)  ! # of available command strings
      parameter (c_ncel= 26)   ! # of available contact elements
      parameter (c_lp1 = 200)  ! # of available history variables definitions
                               !      for vectors CH1 & CH2
      parameter (c_lp3 = 100)  ! # of available history variables definitions
                               !      for vectors CH3
      parameter (c_lmv = 50)   ! # of available variables to store material data
