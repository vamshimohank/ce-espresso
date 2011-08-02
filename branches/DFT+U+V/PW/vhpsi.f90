!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine vhpsi (npwx, npw, nbnd, psip, hpsi)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the Hubbard potential applied to the electronic
  ! of the current k-point, the result is added to hpsi
  !
  USE kinds,     ONLY : DP
  USE parameters,ONLY : lmaxx
  USE becmod,    ONLY : bec_type, calbec, allocate_bec_type, deallocate_bec_type
  USE ldaU,      ONLY : Hubbard_U, Hubbard_alpha, Hubbard_beta, swfcatom, oatwfc, ilm
  USE lsda_mod,  ONLY : nspin, current_spin
  USE scf,       ONLY : v
  USE ions_base, ONLY : nat, ntyp => nsp, ityp
  USE basis,     ONLY : natomwfc
  USE gvect,     ONLY : gstart
  USE control_flags, ONLY : gamma_only
  USE mp_global, ONLY: intra_pool_comm
  USE mp,        ONLY: mp_sum
  !
  implicit none
  !
  integer, intent (in) :: npwx, npw, nbnd
  complex(DP), intent(in) :: psip (npwx, nbnd)
  complex(DP), intent(inout) :: hpsi (npwx, nbnd)
  !
  integer :: ibnd, na, nt, m1, m2, il
  complex(DP) :: temp
  type (bec_type) :: proj
  !
  ! Offset of atomic wavefunctions initialized in setup and stored in oatwfc
  !
  call allocate_bec_type ( natomwfc, nbnd, proj )
  CALL calbec (npw, swfcatom, psip, proj)

  do ibnd = 1, nbnd

    do na = 1, nat  
      nt = ityp (na)  
      do il = 0, lmaxx
        if (oatwfc(na,il) == -1) cycle

        do m1 = 1, 2*il+1
          temp = 0.d0
          if (gamma_only) then
            do m2 = 1, 2*il+1
               temp = temp + v%ns(ilm(il,m1), ilm(il,m2), current_spin, na) * &
                                    proj%r(oatwfc(na,il)+m2, ibnd)
             enddo
             call daxpy (2*npw, temp, swfcatom(1,oatwfc(na,il)+m1), 1, hpsi(1,ibnd), 1)
           else
             do m2 = 1, 2*il+1
               temp = temp + v%ns(ilm(il,m1), ilm(il,m2), current_spin, na) * &
                                    proj%k(oatwfc(na,il)+m2, ibnd)
             enddo
             call zaxpy (npw, temp, swfcatom(1,oatwfc(na,il)+m1), 1, hpsi(1,ibnd), 1)
           endif
         enddo
       enddo
     enddo

  enddo
  call deallocate_bec_type (proj)
  return

end subroutine vhpsi

