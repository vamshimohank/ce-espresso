!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
MODULE el_phon
  USE kinds, ONLY :  DP
  !
  SAVE
  !
  INTEGER :: elph_nbnd_min, elph_nbnd_max
  INTEGER :: el_ph_ngauss, el_ph_nsigma
  REAL(DP) :: el_ph_sigma
  COMPLEX(DP), ALLOCATABLE, TARGET :: &
       el_ph_mat(:,:,:,:)    !  nbnd, nbnd, nks, 3*nat
  !
END MODULE el_phon