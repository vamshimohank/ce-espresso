!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine init_ns
   !-----------------------------------------------------------------------
   !
   ! This routine computes the starting ns (for lda+U calculation) filling
   ! up the d states (the only interested by the on-site potential for the
   ! moment) according to the Hund's rule (valid for the isolated atoms on
   ! which starting potential is built), and to the starting_magnetization:
   ! majority spin levels are populated first, then the remaining electrons
   ! are equally distributed among the minority spin states
   !
   USE kinds,     ONLY : DP
   USE parameters,ONLY : lmaxx
   USE ions_base, ONLY : nat, ityp
   USE lsda_mod,  ONLY : nspin, starting_magnetization
   USE scf,       ONLY : rho
   USE uspp_param,ONLY : upf
   !
   implicit none

   real(DP) :: totoc
   real(DP), external :: hubbard_occ
   integer, external :: set_hubbard_l
   integer :: na, nt, is, m1, majs, mins
   integer :: il, i0, i1, i, lupf
   logical :: magnetic

   rho%ns(:,:,:,:) = 0.d0

   ! loop over atoms
   do na = 1, nat
     nt = ityp (na)

     lupf = set_hubbard_l(upf(nt)%psd)

     ! loop over angular momentum
     do il = 0, lmaxx
       i0 = il*il + 1
       i1 = i0 + 2*il

       ! skip higher angular momentum
       if (il > lupf) cycle

       ! fill the inner shells equally
       if (il < lupf) then
         do i = i0, i1
           if (nspin == 1) then
             rho%ns(i,i,1,na) = 2.d0
           else
             rho%ns(i,i,:,na) = 1.d0
           endif
         enddo
         cycle
       endif

       ! fill the open shell according to Hund's rule
       totoc = Hubbard_occ(upf(nt)%psd)

       magnetic = .false.
       if (nspin == 2) then
         if (starting_magnetization(nt) > 0.d0) then
           magnetic = .true.
           majs = 1
           mins = 2
         elseif (starting_magnetization(nt) < 0.d0) then
           magnetic = .true.
           majs = 2
           mins = 1
         endif
       endif

       ! the atom is spin polarized
       if (magnetic) then

         if (totoc > 2*il+1) then
           do i = i0, i1
             rho%ns(i,i,majs,na) = 1.d0
             rho%ns(i,i,mins,na) = (totoc - real(2*il+1,dp)) / real(2*il+1,dp)
           enddo
         else
           do i = i0, i1
             rho%ns(i,i,majs,na) = totoc / real(2*il+1,dp)
           enddo
         endif
       else ! the atom is not spin polarized
         do i = i0, i1
           rho%ns(i,i,1:2,na) = totoc /  2.d0 / real(2*il+1,dp)
         enddo
       endif  ! if (magnetic)

     enddo ! il
   enddo ! na

end subroutine init_ns
