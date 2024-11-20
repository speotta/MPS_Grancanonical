! Routines need by the Davidson eigensolver (jdqz)
! The only one used are the AMUL routine that performs the multiplication of a sparse matrix with a vector
! and the Guess routine to set the initial guess of the iterative diagonalization.

SUBROUTINE AMUL(n,Psi_Array_In,Psi_Array_Out)

  USE Bose_Hubbard

  IMPLICIT NONE
  INTEGER,                      INTENT(IN)  :: n
  DOUBLE COMPLEX, DIMENSION(n), INTENT(IN)  :: Psi_Array_In
  DOUBLE COMPLEX, DIMENSION(n), INTENT(OUT) :: Psi_Array_Out
  DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: PsiS_In, PsiS_Out
  DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: PsiD_In, PsiD_Out
  
  IF (double_Site_Optimization) THEN
    IF (ldim*rdim*cdim**2 /= n) STOP 'Error in AMUL'
    ALLOCATE(PsiD_In(ldim,rdim,cdim,cdim),PsiD_Out(ldim,rdim,cdim,cdim))
    CALL Array_To_Tensor4(ldim,rdim,cdim,cdim,PsiD_In,n,Psi_Array_In)
    CALL Multiply_Ham_Psi_Two_Site(current_site,ldim,rdim,cdim,PsiD_In,PsiD_Out,Ham,Eff_Ham)
    CALL Tensor4_To_Array(ldim,rdim,cdim,cdim,PsiD_Out,n,Psi_Array_Out)
    DEALLOCATE(PsiD_In,PsiD_Out)
  ELSE
    IF (ldim*rdim*cdim /= n) STOP 'Error in AMUL'
    ALLOCATE(PsiS_In(ldim,rdim,cdim),PsiS_Out(ldim,rdim,cdim))
    CALL Array_To_Tensor3(ldim,rdim,cdim,PsiS_In,n,Psi_Array_In)
    CALL Multiply_Ham_Psi_One_Site(current_site,ldim,rdim,cdim,PsiS_In,PsiS_Out,Ham,Eff_Ham)
    CALL Tensor3_To_Array(ldim,rdim,cdim,PsiS_Out,n,Psi_Array_Out)
    DEALLOCATE(PsiS_In,PsiS_Out)
  END IF

END SUBROUTINE AMUL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BMUL(n,q,r)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE COMPLEX, DIMENSION(n) :: q, r
    
  r=q

END SUBROUTINE BMUL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE PRECON(n,psi)  ! NOT IMPLEMENTED
  !     This subroutine computes $\psi = K \psi$, where $K$ is a matrix which mimics the action of 
  !     $(H - \tau \mathbb{1})^{-1}$. Here H is the Hamiltonian to be diagonalized, $\tau$ is the target 
  !     of the Davidson method, namely the value near which the eigenvalues are sought.
  !     A zeroth order approximation is typically used: $K_{i,j} = \delta_{i,j} (H_{i,i}-\tau)^{-1}$
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE COMPLEX, DIMENSION(n) :: psi

  psi=psi

END SUBROUTINE PRECON

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Guess(n,X)

  USE Bose_Hubbard

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE COMPLEX, DIMENSION(n) :: X

  X = Initial_Guess

END SUBROUTINE Guess

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

