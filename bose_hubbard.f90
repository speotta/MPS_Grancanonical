MODULE Bose_Hubbard
  ! Module with routine specific for the Bose-Hubbard model
  ! and for the static optimization and time evolution
  
  USE Open_DMRG
  
  TYPE(Hamiltonian) Ham						! Hamiltonian fro the static optimization
  TYPE(Effective_Hamiltonian) Eff_Ham				! effective Hamiltonian (MPS contracted on the Hamiltonian)
  DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:) :: Initial_Guess	! initial guess for the optimization
  INTEGER ldim		! link dimension at the left of the current optimization site  
  INTEGER rdim		! link dimension at the right of the current optimization site (double site if double_Site_Optimization is true)
  INTEGER cdim		! dimension of the current optimization site (local Hilbert dimension)
  INTEGER ham_dim	! Effective hamiltonian dimension
  INTEGER current_site	! current site index
  LOGICAL double_Site_Optimization	! flag to signal single or double site optimization
  
  CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Gutzwiller_Translational_Invariant_Guess(u,mu,d,psi)
  ! Produces a traslational invariant Gutzwiller wavefunction for the initial guess

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: d	! dimension of the local Hilbert space 	
  DOUBLE PRECISION, INTENT(IN) :: u	! Hubbard interaction
  DOUBLE PRECISION, INTENT(IN) :: mu	! chemical potential
  DOUBLE PRECISION, DIMENSION(d), INTENT(OUT) :: psi	! local wavefunction that encodes the traslational invariant state
  
  ! working variables
  DOUBLE PRECISION, DIMENSION(d) :: psi_prev, root_n, eigenval, nn
  DOUBLE PRECISION, DIMENSION(d,d) :: ham
  DOUBLE PRECISION :: beta, n
  INTEGER ii, jj

  DO ii = 1, d
    nn(ii) = ii-1.
  ENDDO
  root_n = SQRT(nn)
  psi = 1./SQRT(1.0d0*d)
  DO jj = 1, 1000
    psi_prev = psi
      beta = DOT_PRODUCT(psi,EOSHIFT(root_n*psi,SHIFT = 1))
      ham = 0.
      DO ii = 1, d
	ham(ii,ii) = 0.5*u*(ii-1)*(ii-2) - mu*(ii-1)
      ENDDO
      DO ii = 1, d-1
	ham(ii,ii+1) = -2*beta*root_n(ii+1)
	ham(ii+1,ii) = ham(ii,ii+1)
      ENDDO
      CALL Diagonalize_Real_Symm_Matrix(d,ham,eigenval)
      psi = ham(:,1)
      IF ( ALL(ABS(psi - psi_prev) <= (1.0e-12 + 1.0e-8 * ABS(psi_prev))) ) THEN
	EXIT
      ENDIF
    ENDDO
    n = 0.
    DO ii = 1, d
      n = n + (ii-1)*psi(ii)**2
    ENDDO
    12 FORMAT ('psi = ',20F20.15)
    13 FORMAT ('n   = ',F20.15)
    WRITE(*,12) psi
    WRITE(*,13) n
  !STOP
  CONTAINS
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  SUBROUTINE Diagonalize_Real_Symm_Matrix(d,EigenVec,Eigenval)
    ! Diagonalize a real symmetric matrix stored in EigenVec with dimension d
    ! Store the eigenvectors in EigenVec and the eigenvalues in EigenVal

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: d 	! matrix dimension
    DOUBLE PRECISION, DIMENSION(d,d), INTENT(INOUT) :: EigenVec	! Input: matrix to be diagonalized; Output: eigenvectors
    DOUBLE PRECISION, INTENT(OUT),   DIMENSION(d) :: EigenVal	! Output: eigevalues

    ! working variables
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: Work
    INTEGER :: LWork, Info, ii

    LWork = 3*d-1
    ALLOCATE (Work(LWork))
    CALL DSYEV('V', 'U', d, EigenVec, d, EigenVal, Work, LWork, Info)
    IF (Info /= 0) STOP 'Problems with exact diagonalization'
    DEALLOCATE(Work)
 
  END SUBROUTINE Diagonalize_Real_Symm_Matrix

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
END SUBROUTINE Gutzwiller_Translational_Invariant_Guess

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Build_Number_Operator(d,n)
    ! Build the number operator for a single site
    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: d	! local Hilbert space dimension
    TYPE(One_Site_Operator), INTENT(OUT) :: n	! Output: single site number operator 
    
    ! working variables
    INTEGER ii
    
    n%num_Elements = d
    ALLOCATE(n%element(n%num_Elements))
    DO ii = 1, n%num_Elements
      n%element(ii)%r = ii
      n%element(ii)%c = ii
      n%element(ii)%v = 1.0d0*(ii-1)
    ENDDO

  END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Initialize_MPS_Bose_Hubbard(d,psi,l,MPS)
  ! Initialize a MPS with the translational invariant Gutzwiller wavefunction given in psi

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: l	! lattice length
  INTEGER, INTENT(IN) :: d	! local Hilbert space dimension
  DOUBLE PRECISION, DIMENSION(d), INTENT(IN) :: psi	! traslational invariant Guzwiller wavefunction
  TYPE(Matrix_Product_State), INTENT(OUT) :: MPS	! Output: Matrix Product State
  
  ! working variables
  INTEGER ii, jj
  
  MPS%lattice_length = l
  MPS%local_dim = d
  ALLOCATE(MPS%link_dim(0:MPS%lattice_length))
  MPS%link_dim = 1
  ALLOCATE(MPS%site(MPS%lattice_length))
  DO ii = 1, MPS%lattice_length
    ALLOCATE(MPS%site(ii)%A(MPS%local_dim))
    DO jj = 1, MPS%local_dim
      ALLOCATE(MPS%site(ii)%A(jj)%M(MPS%link_dim(ii-1),MPS%link_dim(ii)))
      MPS%site(ii)%A(jj)%M = psi(jj)
    ENDDO
  ENDDO
  
END SUBROUTINE Initialize_MPS_Bose_Hubbard

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Build_Bose_Hubbard_Hamiltonian(l,d,u,mu,phi,V,Ham)
  ! Build the Bose-Hubbard Hamiltonian
  ! Can be generalized to arbitary Hamiltonians
  ! The Hamiltonian is a collection of l local operators of type One_Site_Operator 
  ! and l-1 nearest-neighbour operators of type Two_Site_Operator
  ! Only the nonzero matrix elements of the local and nearest-neighbour operators are stored

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: l	! lattice length
  INTEGER, INTENT(IN) :: d	! local Hilbert space dimension
  DOUBLE PRECISION, INTENT(IN) :: u	! Hubbard interaction
  DOUBLE PRECISION, INTENT(IN) :: mu	! chemical potential
  DOUBLE PRECISION, INTENT(IN) :: phi	! phase of the hopping term
  DOUBLE PRECISION, DIMENSION(l), INTENT(IN) :: V	! external potential
  TYPE(Hamiltonian), INTENT(OUT) :: Ham			! Output: Hamiltonian

  ! working variables
  DOUBLE PRECISION, DIMENSION(d) :: root_n
  INTEGER ii, jj, kk, counter
  
  DO ii = 1, d
    root_n(ii) = SQRT(ii*1.d0)
  ENDDO
  Ham%lattice_length = l
  Ham%local_dim = d
  ! Local terms of the Hamiltonian
  ALLOCATE(Ham%Local_Terms(Ham%lattice_length))
  DO ii = 1, Ham%lattice_length
    Ham%Local_Terms(ii)%num_Elements = Ham%local_dim
    ALLOCATE(Ham%Local_Terms(ii)%element(Ham%Local_Terms(ii)%num_Elements))
    DO jj = 1, Ham%Local_Terms(ii)%num_Elements
      Ham%Local_Terms(ii)%element(jj)%r = jj
      Ham%Local_Terms(ii)%element(jj)%c = jj
      Ham%Local_Terms(ii)%element(jj)%v = 0.5*u*(jj-1.)*(jj-2.) + (V(ii)-mu)*(jj-1.)
    ENDDO
  ENDDO
  ! Nearest neighbour terms
  ALLOCATE(Ham%Nearest_Neighbour_Terms(Ham%lattice_length-1))
  DO ii = 1, Ham%lattice_length-1
    Ham%Nearest_Neighbour_Terms(ii)%num_Elements = 2*(Ham%local_dim-1)**2
    ALLOCATE(Ham%Nearest_Neighbour_Terms(ii)%element(Ham%Nearest_Neighbour_Terms(ii)%num_Elements))
    counter = 0
    DO jj = 1, Ham%local_dim-1
      DO kk = 1, Ham%local_dim-1
	counter = counter + 1
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%r1 = jj
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%c1 = jj + 1
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%r2 = kk + 1
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%c2 = kk
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%v = -root_n(jj)*root_n(kk)*EXP((0.0,1.0)*phi)
      ENDDO
    ENDDO
    DO jj = 1, Ham%local_dim-1
      DO kk = 1, Ham%local_dim-1
	counter = counter + 1
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%r1 = jj + 1
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%c1 = jj
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%r2 = kk
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%c2 = kk + 1
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%v = -root_n(jj)*root_n(kk)*EXP((0.0,-1.0)*phi)
      ENDDO
    ENDDO
    IF (counter /= Ham%Nearest_Neighbour_Terms(ii)%num_Elements) STOP 'Error in building Bose-Hubbard Hamiltonian'
  ENDDO
  
END SUBROUTINE Build_Bose_Hubbard_Hamiltonian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Build_Bose_Hubbard_Hamiltonian_Time_Evolution(l,d,u,mu,phi,V,Ham)
  ! Build the Bose-Hubbard Hamiltonian for the time evolution 
  ! The local operators are not used, instead they are shared between the nearest neighbour ones
  ! this is convenient since one need to take the exponential of the even and odd parts of the Hamiltonian
  ! H = H_even + H_odd where H_even, H_odd is the sum of mutually commuting two-site operators
  ! The operators comprising H_even (H_odd) are stored in Ham%Nearest_Neighbour_Terms(ii)
  ! where ii is an even (odd) integer
  ! only the nonzero matrix elements of the various operators are stored
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: l		! lattice length
  INTEGER, INTENT(IN) :: d		! local Hilbert space dimension
  DOUBLE PRECISION, INTENT(IN) :: u		! Hubbard interaction
  DOUBLE PRECISION, INTENT(IN) :: mu		! chemical potential
  DOUBLE PRECISION, INTENT(IN) :: phi		! phase of the hopping term
  DOUBLE PRECISION, DIMENSION(l), INTENT(IN) :: V	! external potential
  TYPE(Hamiltonian), INTENT(OUT) :: Ham			! Output: Hamiltonian
  
  ! working variables
  DOUBLE PRECISION, DIMENSION(d) :: root_n, nn
  INTEGER ii, jj, kk, counter
  DOUBLE PRECISION a, b 
  
  DO ii = 1, d
    root_n(ii) = SQRT(ii*1.d0)
  ENDDO
  Ham%lattice_length = l
  Ham%local_dim = d
  ! Nearest neighbour terms
  ALLOCATE(Ham%Nearest_Neighbour_Terms(Ham%lattice_length-1))
  DO ii = 1, Ham%lattice_length-1
    Ham%Nearest_Neighbour_Terms(ii)%num_Elements = Ham%local_dim**2 + 2*(Ham%local_dim-1)**2
    ALLOCATE(Ham%Nearest_Neighbour_Terms(ii)%element(Ham%Nearest_Neighbour_Terms(ii)%num_Elements))
    counter = 0
    IF (ii == 1.) THEN
      a = 1.0
      b = 0.5
    ELSEIF (ii == Ham%lattice_length-1) THEN
      a = 0.5
      b = 1.0
    ELSE
      a = 0.5
      b = 0.5
    ENDIF
    DO jj = 1, Ham%local_dim
      DO kk = 1, Ham%local_dim
	counter = counter + 1
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%r1 = jj
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%c1 = jj
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%r2 = kk
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%c2 = kk
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%v  = a*(0.5*u*(jj-1.)*(jj-2.) + (V(ii)-mu)*(jj-1.)) +  &
							      b*(0.5*u*(kk-1.)*(kk-2.) + (V(ii+1)-mu)*(kk-1.))
      ENDDO
    ENDDO    
    DO jj = 1, Ham%local_dim-1
      DO kk = 1, Ham%local_dim-1
	counter = counter + 1
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%r1 = jj
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%c1 = jj + 1
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%r2 = kk + 1
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%c2 = kk
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%v = -root_n(jj)*root_n(kk)*EXP((0.0,1.0)*phi)
      ENDDO
    ENDDO
    DO jj = 1, Ham%local_dim-1
      DO kk = 1, Ham%local_dim-1
	counter = counter + 1
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%r1 = jj + 1
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%c1 = jj
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%r2 = kk
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%c2 = kk + 1
	Ham%Nearest_Neighbour_Terms(ii)%element(counter)%v = -root_n(jj)*root_n(kk)*EXP((0.0,-1.0)*phi)
      ENDDO
    ENDDO
    IF (counter /= Ham%Nearest_Neighbour_Terms(ii)%num_Elements) STOP 'Error in building Bose-Hubbard Hamiltonian'
  ENDDO
  
END SUBROUTINE Build_Bose_Hubbard_Hamiltonian_Time_Evolution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Optimize_MPS_One_Site(MPS,pos,ground_state_erg)
  ! Optimize a MPS in position pos and return the updated MPS and the ground state energy obtained from the matrix diagonalization
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: pos	! site to be optimized
  TYPE(Matrix_Product_State), INTENT(INOUT) :: MPS	! input-Output: MPS to be optimized
  DOUBLE PRECISION, INTENT(OUT) :: ground_state_erg	! ground state energy with the new MPS
  
  ! working variables
  DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:) :: ground_state_wf
  DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: Temp
  DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: H
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: W
  INTEGER ii
  LOGICAL :: davidson_optimization = .TRUE.	! this flags specifies whether to use the recursive Daivdson algorithm for matrix diagonalization
						! or the full matrix diagonalization. The Davidson is the default, use the other for debugging purposes
  
  IF (pos >= 1 .and. pos <= MPS%lattice_length .and. davidson_optimization) THEN
    ! Update variables in the module (awkward, will be improved in the future)
    current_site = pos
    double_Site_Optimization = .FALSE.
    ldim = MPS%link_dim(current_site-1)
    rdim = MPS%link_dim(current_site)
    cdim = MPS%local_dim
    ham_dim = ldim*rdim*cdim
    !!!!!!!!
    ! Perform the diagonalization, the need to map the MPS to an array stems from the fact that the Davidson routine 
    ! is written in Fortran77 and the interface to the DMRG is not as clean as it should be (may be improved in the future)
    ALLOCATE(Initial_Guess(ham_dim),ground_state_wf(ham_dim),Temp(ldim,rdim,cdim))
    DO ii = 1, cdim
      Temp(:,:,ii) = MPS%site(current_site)%A(ii)%M
    ENDDO
    CALL Tensor3_to_Array(ldim,rdim,cdim,Temp,ham_dim,Initial_Guess)
    CALL Davidson(ground_state_erg,ham_dim,ground_state_wf)
    CALL Array_to_Tensor3(ldim,rdim,cdim,Temp,ham_dim,ground_state_wf)
    DO ii = 1, cdim
      MPS%site(current_site)%A(ii)%M = Temp(:,:,ii)
    ENDDO    
    DEALLOCATE(ground_state_wf,Initial_Guess,Temp)
    !!!!!!!!!!!!!
  ELSEIF (pos >= 1 .and. pos <= MPS%lattice_length) THEN
    ! Update variables in the module (awkward, will be improved in the future)
    current_site = pos
    double_Site_Optimization = .FALSE.
    ldim = MPS%link_dim(current_site-1)
    rdim = MPS%link_dim(current_site)
    cdim = MPS%local_dim
    ham_dim = ldim*rdim*cdim
    !!!!!!!!
    ! Perform the diagonalization using a standard full-matrix diagonalization routine (DiagMatrix)
    ALLOCATE(ground_state_wf(ham_dim),H(ham_dim,ham_dim),Temp(ldim,rdim,cdim),W(ham_dim))
    CALL Build_Hamiltonian_Matrix(ham_dim,H,pos,ldim,rdim,cdim,Ham,Eff_Ham)
    CALL DiagMatrix(ham_dim,H,W)
    ground_state_wf = H(:,1)
    CALL Array_to_Tensor3(ldim,rdim,cdim,Temp,ham_dim,ground_state_wf)
    DO ii = 1, cdim
      MPS%site(current_site)%A(ii)%M = Temp(:,:,ii)
    ENDDO
    ground_state_erg = W(1)
    DEALLOCATE(ground_state_wf,Temp,W)    
  ENDIF

  CONTAINS
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  SUBROUTINE DiagMatrix(n,M,W) ! Uses lapack routines to diagonalize Hermitian matrices

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n	! matrix dimension
  DOUBLE COMPLEX, DIMENSION(n,n), INTENT(INOUT) :: M	! matrix to be diagonalized and output eigenvectors 
  DOUBLE PRECISION, INTENT(OUT), DIMENSION(n) :: W	! eigenvalues
  
  ! working variables
  DOUBLE COMPLEX,   DIMENSION(:),   ALLOCATABLE :: Work
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: RWork
  INTEGER :: LWork, LRWork, Info

  LWork  = 18*n
  LRWork = 3*n-2
  ALLOCATE (Work(LWork),RWork(LRWork))
  CALL ZHEEV('V', 'U', n, M, n, W, Work, LWork, RWork, Info)
  IF (Info /= 0) STOP 'Problems with exact diagonalization'
  DEALLOCATE(Work,RWork)

  END SUBROUTINE DiagMatrix
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
END SUBROUTINE Optimize_MPS_One_Site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Optimize_MPS_Two_Site(MPS,pos,ground_state_erg,m_Min,m_Max,eps,LeftRight)
  ! Two-site optimization of the MPS on position pos
  ! First finds the ground state of the effective Hamiltonian then truncate using the SVD
  
  IMPLICIT NONE
  CHARACTER, INTENT(IN) :: LeftRight	! Flag indicating whether to absorb the singular values in the left MPS matrices (right gauge) or on the right MPS matrices (left gauge)
  INTEGER, INTENT(IN) :: pos		! site to be optimized
  INTEGER, INTENT(IN) :: m_Min		! minimum value of the link dimension
  INTEGER, INTENT(IN) :: m_Max		! maximum value of the link dimension
  DOUBLE PRECISION, INTENT(IN) :: eps	! truncation error
  TYPE(Matrix_Product_State), INTENT(INOUT) :: MPS	! Matrix Product State to be optimized
  DOUBLE PRECISION, INTENT(OUT) :: ground_state_erg	! ground state energy obtained from the diagonalization of the effective Hamiltonian

  ! working variables
  DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:) :: ground_state_wf
  DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: Temp
  DOUBLE COMPLEX, DIMENSION(:,:), POINTER :: B
  DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE :: U, V
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: S, P
  INTEGER ii, jj, dimMin
  DOUBLE PRECISION ErrorSum
  
  IF (pos >= 1 .and. pos < MPS%lattice_length) THEN
    ! Update variables in the module (awkward, will be improved in the future)
    current_site = pos
    double_Site_Optimization = .TRUE.
    ldim = MPS%link_dim(current_site-1)
    rdim = MPS%link_dim(current_site+1)
    cdim = MPS%local_dim
    ham_dim = ldim*rdim*cdim**2
    !!!!!!!!!
    ! Perform the diagonalization, the need to map the MPS to an array stems from the fact that the Davidson routine 
    ! is written in Fortran77 and the interface to the DMRG is not as clean as it should be (may be improved in the future)
    ALLOCATE(Initial_Guess(ham_dim),ground_state_wf(ham_dim),Temp(ldim,rdim,cdim,cdim))
    DO ii = 1, cdim
      DO jj = 1, cdim
	Temp(:,:,ii,jj) = MATMUL(MPS%site(current_site)%A(ii)%M,MPS%site(current_site+1)%A(jj)%M)
      ENDDO    
    ENDDO
    CALL Tensor4_to_Array(ldim,rdim,cdim,cdim,Temp,ham_dim,Initial_Guess)
    CALL Davidson(ground_state_erg,ham_dim,ground_state_wf)
    CALL Array_to_Tensor4(ldim,rdim,cdim,cdim,Temp,ham_dim,ground_state_wf)
    ! Perform the SVD
    ALLOCATE(B(ldim*cdim,rdim*cdim))
    B = 0.
    DO ii = 1, cdim
      DO jj = 1, cdim 
	B(1+(ii-1)*ldim:ii*ldim,1+(jj-1)*rdim:jj*rdim) = Temp(:,:,ii,jj)
      ENDDO
    ENDDO
    dimMin = MIN(ldim,rdim)*cdim
    ALLOCATE(U(ldim*cdim,dimMin),V(dimMin,rdim*cdim),S(dimMin))
    CALL SVD(ldim*cdim,dimMin,rdim*cdim,B,U,S,V)
    ! Truncate
    ALLOCATE(P(dimMin))
    P = S**2/SUM(S**2)
    ErrorSum = 0.
    IF (dimMin > m_Min) THEN
      loop1: DO ii = dimMin, 1, -1
	  IF ( ii <= m_Max .AND. ((ErrorSum + P(ii)) > eps) .OR. ii <= m_Min) THEN
	    MPS%link_dim(pos) = ii
	    EXIT loop1
	  ENDIF
	  ErrorSum = ErrorSum + P(ii)
      ENDDO loop1
    ENDIF
    ! Absorb singular values on the left or right site according to LeftRight
    IF (LeftRight == 'L') THEN ! Absorb on the left site (Right gauge)
      DO ii = 1, MPS%link_dim(pos)
	U(:,ii) = U(:,ii)*S(ii) 
      ENDDO
    ELSEIF (LeftRight == 'R') THEN ! Absorb on the right site (Left gauge)
      DO ii = 1, MPS%link_dim(pos)
	V(ii,:) = S(ii)*V(ii,:) 
      ENDDO  
    ENDIF
    ! Update the MPS with the new approximate ground state
    DO ii = 1, cdim
      DEALLOCATE(MPS%site(pos)%A(ii)%M)
      DEALLOCATE(MPS%site(pos+1)%A(ii)%M)
      ALLOCATE(MPS%site(pos)%A(ii)%M(ldim,MPS%link_dim(pos)))
      ALLOCATE(MPS%site(pos+1)%A(ii)%M(MPS%link_dim(pos),rdim))
      MPS%site(pos)%A(ii)%M = U(1+(ii-1)*ldim:ii*ldim,1:MPS%link_dim(pos))
      MPS%site(pos+1)%A(ii)%M = V(1:MPS%link_dim(pos),1+(ii-1)*rdim:ii*rdim)
    ENDDO
    DEALLOCATE(B,U,S,V,ground_state_wf,Initial_Guess,Temp,P)
  ENDIF
  
  CONTAINS
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  SUBROUTINE SVD(dimRow,dimMin,dimCol,M,U,S,V)
  ! Perform the Singular Value Decompostion
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: dimRow		! number of rows
  INTEGER, INTENT(IN) :: dimCol		! number of columns
  INTEGER, INTENT(IN) :: dimMin		! minimum between rows and column dimension
  DOUBLE COMPLEX, DIMENSION(dimRow,dimCol), INTENT(IN)  :: M	! Input: matrix to be decomposed
  DOUBLE COMPLEX, DIMENSION(dimRow,dimMin), INTENT(OUT) :: U	! Output: Left unitray matrix
  DOUBLE COMPLEX, DIMENSION(dimMin,dimCol), INTENT(OUT) :: V	! Output: Right unitary matrix
  DOUBLE PRECISION, DIMENSION(dimMin), INTENT(OUT) :: S		! Output: singular values
  
  ! working variables
  DOUBLE COMPLEX,   DIMENSION(:), ALLOCATABLE :: Work
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RWork
  INTEGER  :: LWork, LRWork, Info
  
  LWork  = 24 * (2*dimMin + Max(dimRow,dimCol))
  LRWork = 5  * dimMin
  ALLOCATE(Work(LWork),RWork(LRWork))
  CALL ZGESVD('S','S',dimRow,dimCol,M,dimRow,S,U,dimRow,V,dimMin,Work,LWork,RWork,Info)
  IF (Info /= 0) STOP 'Problems in SVD calculation'
  
  END SUBROUTINE SVD
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
END SUBROUTINE Optimize_MPS_Two_Site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
SUBROUTINE MPS_Jacobi_Davidson_Optimization_Single_Site(MPS,sweeps_Single_Site)
  ! Perform several sweeps of single site optimization on a given MPS
  ! the Hamiltonian is specified by Ham and related variables in this module
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: sweeps_Single_Site	! number of sweeps
  TYPE(Matrix_Product_State), INTENT(INOUT) :: MPS		! Matrix Product State
  
  ! working variables
  DOUBLE PRECISION ground_state_erg
  INTEGER sweep, ii
  REAL time1, time2
  
  IF (sweeps_Single_Site > 0) THEN
  ! One sweeps corresponds to optimizing from site 1 to site l and back
  DO sweep = 1, sweeps_Single_Site
    ! Sweep from site 1 to site l-1
    DO ii = 1, MPS%lattice_length-1
      CALL CPU_time(time1)
      ! Optimize site ii
      CALL Optimize_MPS_One_Site(MPS,ii,ground_state_erg)
      CALL CPU_time(time2)
      print *, ii, ham_dim, ground_state_erg, time2-time1
      WRITE(104,*) ground_state_erg
      ! Update the gauge to optimize in the next site
      CALL MPS_Singular_Value_Decomposition(MPS,ii,'R')
      ! Update the Effective Hamiltonian with the new MPS
      CALL Update_Effective_Hamiltonian(MPS,Ham,Eff_Ham,ii+1,'L')
    ENDDO
    ! Sweep from site l to site 2
    DO ii = MPS%lattice_length, 2, -1
      CALL CPU_time(time1)
      ! Optimize site ii
      CALL Optimize_MPS_One_Site(MPS,ii,ground_state_erg)
      CALL CPU_time(time2)
      print *, ii, ham_dim, ground_state_erg, time2-time1
      WRITE(104,*) ground_state_erg
      ! Update the gauge to optimize in the next site
      CALL MPS_Singular_Value_Decomposition(MPS,ii,'L')
      ! Update the Effective Hamiltonian with the new MPS      
      CALL Update_Effective_Hamiltonian(MPS,Ham,Eff_Ham,ii-1,'R')
    ENDDO
  ENDDO
  ! Final optimization on site 1
  CALL CPU_time(time1)
  CALL Optimize_MPS_One_Site(MPS,1,ground_state_erg)
  CALL CPU_time(time2)
  print *, 1, ham_dim, ground_state_erg, time2-time1
  WRITE(104,*) ground_state_erg
  ENDIF
  
END SUBROUTINE MPS_Jacobi_Davidson_Optimization_Single_Site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MPS_Jacobi_Davidson_Optimization_White_Correction (MPS,sweeps_White_Correction,m_Min,m_Max,eps)
  ! Perform several sweeps of single site optimization on a given MPS with 
  ! the Hamiltonian specified by Ham and related variables in this module
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: sweeps_White_Correction	! number of sweeps
  INTEGER, INTENT(IN) :: m_Min		! minimum value of the link dimension
  INTEGER, INTENT(IN) :: m_Max		! maximum value of the link dimension
  DOUBLE PRECISION, INTENT(IN) :: eps	! truncation error
  TYPE(Matrix_Product_State), INTENT(INOUT) :: MPS 	! Matrix Product State to be optimized
  
  ! working variables
  DOUBLE PRECISION ground_state_erg
  INTEGER sweep, ii
  REAL time1, time2
  DOUBLE PRECISION, PARAMETER :: White_Correction_Magnitude = 0.1  ! parameter controlling the White correction
 
  IF (sweeps_White_Correction > 0) THEN
  ! One sweeps corresponds to optimizing from site 1 to site l and back
  DO sweep = 1, sweeps_White_Correction
    ! Sweep from site 1 to site l-1
    DO ii = 1, MPS%lattice_length-1
      CALL CPU_time(time1)
      ! Optimize site ii
      CALL Optimize_MPS_One_Site(MPS,ii,ground_state_erg)
      CALL CPU_time(time2)
      print *, ii, ham_dim, ground_state_erg, time2-time1
      WRITE(104,*) ground_state_erg
      ! Add to the MPS the nearest neighbour term relative to site ii and ii+1 
      ! multiplied by White_Correction_Magnitude and truncate
      CALL MPS_Multiply_and_Truncate(MPS,Ham%Nearest_Neighbour_Terms(ii), &
		ii,m_Min,m_Max,eps,'R',coeff=White_Correction_Magnitude)
      ! Update the effective Hamiltonian
      CALL Update_Effective_Hamiltonian(MPS,Ham,Eff_Ham,ii+1,'L')
    ENDDO
    ! Sweep from site l to site 2
    DO ii = MPS%lattice_length, 2, -1
      CALL CPU_time(time1)
      ! Optimize site ii
      CALL Optimize_MPS_One_Site(MPS,ii,ground_state_erg)
      CALL CPU_time(time2)
      print *, ii, ham_dim, ground_state_erg, time2-time1
      WRITE(104,*) ground_state_erg
      ! Add to the MPS the nearest neighbour term relative to site ii and ii+1 
      ! multiplied by White_Correction_Magnitude and truncate
      CALL MPS_Multiply_and_Truncate(MPS,Ham%Nearest_Neighbour_Terms(ii-1), &
		ii-1,m_Min,m_Max,eps,'L',coeff=White_Correction_Magnitude)
      ! Update the effective Hamiltonian
      CALL Update_Effective_Hamiltonian(MPS,Ham,Eff_Ham,ii-1,'R')
    ENDDO
  ENDDO
  ! Final optimization on site 1
  CALL CPU_time(time1)
  CALL Optimize_MPS_One_Site(MPS,1,ground_state_erg)
  CALL CPU_time(time2)
  print *, 1, ham_dim, ground_state_erg, time2-time1
  WRITE(104,*) ground_state_erg
  ENDIF
  
END SUBROUTINE MPS_Jacobi_Davidson_Optimization_White_Correction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MPS_Jacobi_Davidson_Optimization_Double_Site(MPS,sweeps_Double_Site,m_Min,m_Max,eps)
  ! Perform several sweeps of double site optimization on a given MPS with 
  ! the Hamiltonian specified by Ham and related variables in this module 

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: sweeps_Double_site	! number of sweeps
  INTEGER, INTENT(IN) :: m_Min		! minimum value of the link dimension
  INTEGER, INTENT(IN) :: m_Max		! maximum value of the link dimension
  DOUBLE PRECISION, INTENT(IN) :: eps	! truncation error
  TYPE(Matrix_Product_State), INTENT(INOUT) :: MPS 	! Matrix Product State to be optimized

  ! working variables
  DOUBLE PRECISION ground_state_erg
  INTEGER sweep, ii
  REAL time1, time2

  IF (sweeps_Double_Site > 0) THEN
  ! One sweeps corresponds to optimizing from site 1 to site l and back
  DO sweep = 1, sweeps_Double_Site
    ! Sweep from site 1 to site l-1    
    DO ii = 1, MPS%lattice_length-2
      CALL CPU_time(time1)
      ! Oprtimize site ii and ii+1 and perform the SVD
      ! the singular values are absorbed on the right site (ii+1)
      CALL Optimize_MPS_Two_Site(MPS,ii,ground_state_erg,m_Min,m_Max,eps,'R')
      CALL CPU_time(time2)
      print *, ii, ham_dim, ground_state_erg, time2-time1
      WRITE(104,*) ground_state_erg
      ! Update the effective Hamiltonian
      CALL Update_Effective_Hamiltonian(MPS,Ham,Eff_Ham,ii+1,'L')
    ENDDO
    ! Sweep from site l to site 2
    DO ii = MPS%lattice_length-1, 2, -1
      CALL CPU_time(time1)
      ! Oprtimize site ii and ii+1 and perform the SVD
      ! the singular values are absorbed on the left site (ii)
      CALL Optimize_MPS_Two_Site(MPS,ii,ground_state_erg,m_Min,m_Max,eps,'L')
      CALL CPU_time(time2)
      print *, ii, ham_dim, ground_state_erg, time2-time1
      WRITE(104,*) ground_state_erg
      ! Update the effective Hamiltonian
      CALL Update_Effective_Hamiltonian(MPS,Ham,Eff_Ham,ii,'R')
    ENDDO
  ENDDO
  CALL CPU_time(time1)
  ! Optimize site ii and ii+1
  CALL Optimize_MPS_Two_Site(MPS,1,ground_state_erg,m_Min,m_Max,eps,'L')
  CALL CPU_time(time2)
  print *, 1, ham_dim, ground_state_erg, time2-time1
  WRITE(104,*) ground_state_erg
  ! Update the effective Hamiltonian
  CALL Update_Effective_Hamiltonian(MPS,Ham,Eff_Ham,1,'R')
  ENDIF
  
END SUBROUTINE MPS_Jacobi_Davidson_Optimization_Double_Site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE TimeEvolution(MPS,time_Ham,Imaginary,Trotter_Order,step,time_steps,obs_Steps,m_Min,m_Max,eps)
  ! Subroutine for performing real and imaginary time-evolution

  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: Imaginary	! Flag to specify either the real or the imaginary time evolution
  INTEGER, INTENT(IN) :: time_steps	! number time steps
  INTEGER, INTENT(IN) :: Trotter_Order	! order of the Trotter expansion
  INTEGER, INTENT(IN) :: obs_steps	! number of steps between the measure of the observables
  INTEGER, INTENT(IN) :: m_Min		! minimum link dimension
  INTEGER, INTENT(IN) :: m_Max		! maximum link dimension 
  DOUBLE PRECISION, INTENT(IN) :: eps	! truncation error
  DOUBLE PRECISION, INTENT(IN) :: step	! time step
  TYPE(Matrix_Product_State), INTENT(INOUT) :: MPS	! Matrix Product State to be evolved
  TYPE(Hamiltonian), INTENT(IN) :: time_Ham		! Hamiltonian used for the time-evolution

  ! working variables
  INTEGER NumberOfTerms, ii, jj, kk
  CHARACTER LeftRight
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: expCoeff
  TYPE(Two_Site_Operator), DIMENSION(:,:), ALLOCATABLE :: expHam
  TYPE(Block_Info_Type) :: binf
  
  DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE :: out
  TYPE(One_Site_Operator) :: n
  
  INTEGER r1, r2, c1, c2
  DOUBLE PRECISION N1, N2, norm1, norm2
  
  21 FORMAT(1000ES10.3)
  22 FORMAT('t = ',F16.5)
  
  SELECT CASE(Trotter_Order)
  CASE(2)
     ! Secondo order expansion
     NumberOfTerms = 3
     IF (ALLOCATED(expCoeff)) DEALLOCATE(expCoeff)
     ALLOCATE(expCoeff(NumberOfTerms))
     expCoeff = [0.5d0, 1.0d0, 0.5d0]
  CASE(4)
     ! Fourth order expansion
     NumberOfTerms = 7
     IF (ALLOCATED(expCoeff)) DEALLOCATE(expCoeff)
     ALLOCATE(expCoeff(NumberOfTerms))
     expCoeff = [0.67560359597982889d0,1.3512071919596578d0,-0.17560359597982883d0,  &
          -1.7024143839193153d0,-0.17560359597982883d0,1.3512071919596578d0,0.67560359597982889d0]
  CASE(6)
     ! Sixth order expansion
     NumberOfTerms = 15
     IF (ALLOCATED(expCoeff)) DEALLOCATE(expCoeff)
     ALLOCATE(expCoeff(NumberOfTerms))
     expCoeff = [0.39225680523877998, 0.78451361047755996, 0.51004341191845848, 0.23557321335935699, &
          -0.47105338540975655, -1.1776799841788701, 0.068753168252518093, 1.3151863206839063,  &
          0.068753168252518093, -1.1776799841788701, -0.47105338540975655, 0.23557321335935699, &
          0.51004341191845848, 0.78451361047755996, 0.39225680523877998]
  CASE(8)
     ! Eight order expansion
     NumberOfTerms = 31
     IF (ALLOCATED(expCoeff)) DEALLOCATE(expCoeff)
     ALLOCATE(expCoeff(NumberOfTerms))
     expCoeff = [0.31451532510521651, 0.62903065021043303, 0.99919005718957155, 1.3693494641687101, &
          0.15238115813844, -1.0645871478918301, 0.29938547587066, 1.6633580996331501, &
          -0.0078055914816249627, -1.6789692825964, -1.6192186604054351, -1.5594680382144701, &
          -0.62383861289802156, 0.31179081241842699, 0.98539084848119352, 1.65899088454396, & 
          0.98539084848119352, 0.31179081241842699, -0.62383861289802156, -1.5594680382144701, & 
          -1.6192186604054351, -1.6789692825964, -0.0078055914816249627, 1.6633580996331501, & 
          0.29938547587066, -1.0645871478918301, 0.15238115813844, 1.3693494641687101, &
          0.99919005718957155, 0.62903065021043303, 0.31451532510521651]
  END SELECT
  
  ! Finds the blocks with nonzero matrix elements in the exponentiated nearest-neighbour interactions
  ! useful to exponentiate only the block and not the whole operators
  ! mainly useful when the Hamiltonian is a fucntion of time and the exponentials need to be calculated at each time step
  CALL Find_Blocks(MPS%local_dim,time_Ham%Nearest_Neighbour_Terms(1),binf)
  ! Calculate the exponentials of the nearest neighbor terms using the info on the blocks obtained before 
  ALLOCATE (expHam(MPS%lattice_length-1,NumberOfTerms))
  DO jj = 1, NumberOfTerms
    DO ii = 1, MPS%lattice_length-1
      IF (Imaginary) THEN
	CALL Exponentiate(time_Ham%Nearest_Neighbour_Terms(ii),binf,(-1.0d0,0.0d0)*step*expCoeff(jj),expHam(ii,jj))
      ELSE
	CALL Exponentiate(time_Ham%Nearest_Neighbour_Terms(ii),binf,(0.0d0,-1.0d0)*step*expCoeff(jj),expHam(ii,jj))
      ENDIF
    ENDDO
  ENDDO

  ! Build the number operator
  call Build_Number_Operator(MPS%local_dim,n)
  
  ALLOCATE(out(MPS%lattice_length))

  ! Normalize the MPS in the right gauge
  CALL Normalize_MPS(MPS,'R')
  
  DO ii = 1, time_steps

    LeftRight = 'R' ! Direction of sweeping of the chain, 'R' from left to right, 'L' from right to left
    DO jj = 1, NumberOfTerms
      IF (LeftRight == 'R') THEN
	! sweep from left to right, H_odd = sum_i h_{i,i+1} with i odd
	DO kk = 1, MPS%lattice_length-1, 2
	  ! Multiply the MPS by the exponential of nearest-neighbor term expHam(kk,jj) and perform SVD
	  ! the truncation can have the particle number correction
	  CALL MPS_Multiply_and_Truncate_with_Particle_Number_Correction(MPS,expHam(kk,jj),kk,m_Min,m_Max,eps,LeftRight,.TRUE.)
	  ! Update the gauge for the next multiplication and truncation
	  CALL MPS_Singular_Value_Decomposition(MPS,kk+1,LeftRight)
	ENDDO
	LeftRight = 'L'
      ELSEIF (LeftRight == 'L') THEN
	! Update the gauge on the final site to prepare for the sweep from right to left
	CALL MPS_Singular_Value_Decomposition(MPS,MPS%lattice_length,'L')
	! sweep from right to left, H_even = sum_i h_{i,i+1} with i even 
	DO kk = MPS%lattice_length-2, 2, -2
	  ! Multiply the MPS by the exponential of nearest-neighbor term expHam(kk,jj) and perform SVD
	  ! the truncation can have the particle number correction
	  CALL MPS_Multiply_and_Truncate_with_Particle_Number_Correction(MPS,expHam(kk,jj),kk,m_Min,m_Max,eps,LeftRight,.TRUE.)
	  ! Update the gauge for the next multiplication and truncation
	  CALL MPS_Singular_Value_Decomposition(MPS,kk,LeftRight)
	ENDDO
	LeftRight = 'R'
      ENDIF
    ENDDO
    ! Normalize the MPS in the right gauge
    CALL Normalize_MPS(MPS,'R')

    ! Measure the observables every 'obs_Steps' steps
    ! only density implemented
    IF (MOD(ii,obs_Steps) == 0) THEN
      print 22, step*ii    
      
      CALL Measure_Local_Observable(MPS,n,out)
      OPEN(unit=103, file='./density.dat', status = 'OLD', position = 'APPEND')
      WRITE(103,21) REAL(out)
      CLOSE(unit = 103)
    ENDIF
  ENDDO
  
!  CONTAINS
  
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
!   SUBROUTINE Test(dd,cd,AA,BB,mm,NN)
!   
!   IMPLICIT NONE
!   INTEGER, INTENT(IN) :: dd, cd
!   TYPE(Matrix), DIMENSION(dd), INTENT(IN) :: AA
!   TYPE(Matrix), DIMENSION(dd), INTENT(IN) :: BB
!   DOUBLE PRECISION, INTENT(OUT) :: mm, NN
!   DOUBLE COMPLEX, DIMENSION(cd,cd) :: Temp1, Temp2, Temp1_N, Temp2_N
!   DOUBLE COMPLEX Normalization
!   INTEGER ii
!   
!   Temp1 = 0.d0
!   Temp2 = 0.d0
!   Temp1_N = 0.d0
!   Temp2_N = 0.d0
!   DO ii = 1, dd
!     Temp1 = Temp1 + MATMUL( TRANSPOSE(CONJG(AA(ii)%M)), AA(ii)%M )
!     Temp2 = Temp2 + MATMUL( BB(ii)%M, TRANSPOSE(CONJG(BB(ii)%M)) )
!     Temp1_N = Temp1_N + (ii-1.0d0)*MATMUL( TRANSPOSE(CONJG(AA(ii)%M)), AA(ii)%M )
!     Temp2_N = Temp2_N + (ii-1.0d0)*MATMUL( BB(ii)%M, TRANSPOSE(CONJG(BB(ii)%M)) )
!   ENDDO
!   mm = Trace(cd,Temp1,Temp2)
!   NN = Trace(cd,Temp1_N,Temp2) + Trace(cd,Temp1,Temp2_N)
!   
!   END SUBROUTINE Test
!   
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
!   PURE FUNCTION Trace(dd,A,B)
!     
!     IMPLICIT NONE
!     INTEGER, INTENT(IN) :: dd
!     DOUBLE COMPLEX, DIMENSION(dd,dd), INTENT(IN) :: A, B
!     DOUBLE COMPLEX Trace
!     INTEGER ii, jj
!     
!     Trace = 0.d0
!     DO ii = 1, dd
!       DO jj = 1, dd
! 	Trace = Trace + A(ii,jj)*B(jj,ii)
!       ENDDO
!     ENDDO
!   
!   END FUNCTION Trace
!   
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
END SUBROUTINE TimeEvolution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE Bose_Hubbard
