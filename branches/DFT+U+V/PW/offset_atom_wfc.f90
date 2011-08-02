!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE offset_atom_wfc( nat, offset )
  !----------------------------------------------------------------------------
  !
  ! For each Hubbard atom, compute the index of the projector in the
  ! list of atomic wavefunctions
  !
  USE uspp_param,       ONLY : upf
  USE parameters,       ONLY : lmaxx
  USE noncollin_module, ONLY : noncolin
  USE ions_base,        ONLY : ityp
  USE basis,            ONLY : natomwfc
  USE mp_global,        ONLY : stdout
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: nat
  INTEGER, INTENT(OUT) :: offset(nat,0:lmaxx)
  INTEGER  :: counter, na, nt, n, il

  counter = 0       ! counting from zero
  offset(:,:) = -1

  write(stdout,*)
  write(stdout,'(5X,''Generating atomic orbitals: -------------------------------------'')')

  do na = 1, nat
     nt = ityp(na)

     do n = 1, upf(nt)%nwfc
        if ( upf(nt)%oc(n) >= 0.d0 ) then
           if (noncolin) then
              call errore('Hubbard_offset_atom_wfc', 'LDA+U non-collinear not yer implemented', 1)
              !!if (upf(nt)%has_so) then
              !!   counter = counter + 2*upf(nt)%lchi(n)
              !!   if (abs(upf(nt)%jchi(n)-upf(nt)%lchi(n) - 0.5d0) < 1.d-6) counter = counter + 2
              !!else
              !!   counter = counter + 2 * ( 2 * upf(nt)%lchi(n) + 1 )
              !!end if
           else
              il = upf(nt)%lchi(n)
              if (offset(na,il) /= -1) &
                write(stdout,'(5X,''warning: atom '',I4,'', l='',I1,'' has two atomic orbitals'')') na, il
              offset(na,il) = counter
              !!write(stdout,'(''na='',I3,4X,''il='',I1,6X,''counter+1='',I4)') na, il, counter+1
              counter = counter + 2*il + 1
           endif
        endif
     enddo
     write(stdout,'(5X,''atom '',I4,4X,''spdf offset:'',4(I4))') na, (offset(na,il), il=0,lmaxx)
  enddo

  write(stdout,'(5X,''-----------------------------------------------------------------'')')
  write(stdout,*)

  if (counter /= natomwfc) call errore('offset_atom_wfc', 'Serious problem: counter /= natomwfc', 1)

END SUBROUTINE offset_atom_wfc

