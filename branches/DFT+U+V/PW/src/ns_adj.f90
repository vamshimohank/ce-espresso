!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine ns_adj
!-----------------------------------------------------------------------
! This routine tries to suggest to the code the right atomic orbital to 
! localize the charge on.
!
   USE kinds,     ONLY : DP
   USE parameters,ONLY : lmaxx
   USE ions_base, ONLY : nat, ntyp => nsp, ityp
   USE ldaU,      ONLY : Hubbard_U, starting_ns
   USE scf,       ONLY : rho
   USE lsda_mod,  ONLY : nspin
   USE io_global, ONLY : stdout
 
   implicit none
   !
   integer, parameter :: ldim = (lmaxx+1)**2
   integer :: na,nt,is,m1,m2,majs,mins,adjs,mol(ldim),nel,i,j,l,index(ldim),il
   real(DP) :: totoc, delta,lambda(ldim) 
   complex(DP) :: vet(ldim,ldim), f(ldim,ldim), temp
   logical :: adjust
 
   if (ALL(starting_ns == -1.d0)) return
   write (stdout,*) "Modify starting ns matrices according to input values "
 
   do na = 1,nat
     nt = ityp(na)
     do il = 0, lmaxx
       if (Hubbard_U(nt,il,il) == 0.d0) cycle
         do is = 1,nspin
           do m1 = 1, 2*il+1
             do m2 = 1, 2*il+1
               f(m1,m2) = rho%ns(il*il+m1,il*il+m2,is,na)
             end do
          end do
          call cdiagh(2*il+1, f, ldim, lambda, vet)
          do i = 1, 2*il+1
            if (starting_ns(il*il+i,is,nt) >= 0.d0) lambda(i) = starting_ns(il*il+i,is,nt)
          end do
          do m1 = 1, 2*il+1
            do m2 = m1, 2*il+1
              temp = 0.d0
              do i = 1, 2*il+1
                 temp = temp + CONJG(vet(m1,i))*lambda(i)*vet(m2,i)     
              end do
              rho%ns(il*il+m1,il*il+m2,is,na) =  DBLE(temp)
              rho%ns(il*il+m2,il*il+m1,is,na) = rho%ns(il*il+m1,il*il+m2,is,na)
            end do
          end do
       end do
     end do
   end do ! on na
   CALL write_ns
   return
end subroutine ns_adj

