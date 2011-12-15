!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine write_ns
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : lmaxx
  USE constants,  ONLY : rytoev
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE lsda_mod,   ONLY : nspin
  USE io_global,  ONLY : stdout
  USE scf,        ONLY : rho
  USE ldaU,       ONLY : Hubbard_U, Hubbard_alpha, Hubbard_beta, oatwfc, ilm
  !
  implicit none
  !
  real(dp) :: tr(0:lmaxx,2)
  real(dp) :: nsum, msum
  integer :: is, il, i0, i1, na, nt, m1
  character(len=2) :: spin_label

  write(stdout,*)
  write(stdout,'(5X,''Occupation matrix for Hubbard: -------------------------------------------'')')
  nsum = 0.d0
  msum = 0.d0

  ! loop over atoms and angular momentum
  do na = 1, nat
    nt = ityp (na)

    ! loop over spin
    do is = 1, nspin
      tr(0:lmaxx,is) = 0.d0

      do il = 0, lmaxx
        ! find offset of atomic wfcs
        if (oatwfc(na,il) == -1) cycle

        do m1 = 1, 2*il+1
          tr(il,is) = tr(il,is) + rho%ns(ilm(il,m1),ilm(il,m1),is,na)
        enddo
        if (nspin == 1) tr(il,1) = tr(il,1) * 2.d0
      enddo ! il

      spin_label = '..'
      if (nspin == 2 .and. is == 1) spin_label = 'up'
      if (nspin == 2 .and. is == 2) spin_label = 'dw'
      write(stdout,'(5X,''atom'',I4,4X,A2,4X,''  s:'',F10.4,''  p:'',F10.4,''  d:'',F10.4,''  f:'',F10.4)') &
            na, spin_label, tr(0:lmaxx,is)
      nsum = nsum + sum(tr(0:lmaxx,is))
      if (nspin == 2 .and. is == 1) msum = msum + sum(tr(0:lmaxx,is))
      if (nspin == 2 .and. is == 2) msum = msum - sum(tr(0:lmaxx,is))
    enddo ! is

  enddo ! na
  write(stdout,'(5X,''Sum of occupations = '',F10.4)') nsum
  if (nspin == 2) write(stdout,'(5X,''Magnetization      = '', F10.4)') msum
  write(stdout,'(5X,''--------------------------------------------------------------------------'')')
  write(stdout,*)

end subroutine write_ns

