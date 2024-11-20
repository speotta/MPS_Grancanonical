MODULE Open_DMRG
  ! Essential data types and subroutines for DMRG and TDMRG
  
  IMPLICIT NONE
  
  ! A matrix with unspecified dimensions
  TYPE Matrix
    DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: M
  END TYPE Matrix
  
  ! Matrix element of a one-site operator
  TYPE One_Site_Operator_Matrix_Element
    INTEGER r, c
    DOUBLE COMPLEX v
  END TYPE One_Site_Operator_Matrix_Element
  
  ! One-site operator: a collections of one-site operator matrix elements
  TYPE One_Site_Operator
    INTEGER num_Elements
    TYPE(One_Site_Operator_Matrix_Element), ALLOCATABLE, DIMENSION(:) :: element
  END TYPE One_Site_Operator
  
  ! Matrix elements of a two-site operator
  TYPE Two_Site_Operator_Matrix_Element
    INTEGER r1, c1, r2, c2
    DOUBLE COMPLEX v
  END TYPE Two_Site_Operator_Matrix_Element
  
  ! Two-site operator: a collection of two-site operator matrix elements
  TYPE Two_Site_Operator
    INTEGER num_Elements
    TYPE(Two_Site_Operator_Matrix_Element), ALLOCATABLE, DIMENSION(:) :: element
  END TYPE Two_Site_Operator
  
  ! Matrix element of a matrix-valued operator, i.e. a matrix
  TYPE Matrix_Valued_Operator_Element
    INTEGER r, c
    DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: M
  END TYPE Matrix_Valued_Operator_Element
  
  ! Matrix-valued operator, i.e. a collecton of matrix-valued operator matrix elements
  TYPE Matrix_Valued_Operator
    INTEGER num_Elements
    TYPE(Matrix_Valued_operator_Element), ALLOCATABLE, DIMENSION(:) :: element
  END TYPE Matrix_Valued_Operator

  ! A collection of d (d is the unspecified local dimension of the Hilbert space) matrices relative to a single site
  ! The building block of a MPS
  TYPE Matrix_Product_State_Site
    TYPE(Matrix), ALLOCATABLE, DIMENSION(:) :: A 
  END TYPE Matrix_Product_State_Site

  ! Matrix product state, a collection of d (d is the unspecified local dimension of the Hilbert space)
  ! matrices for each site (stored in array 'site' of type Matrix_Product_State_Site)
  TYPE Matrix_Product_State
    INTEGER lattice_length 	! lattice length
    INTEGER local_dim		! local dimension of the Hilbert space (also called 'd')
    INTEGER, ALLOCATABLE, DIMENSION(:) :: link_dim	! array storing the link dimensions, the dimension of the MPS matrices
    TYPE(Matrix_Product_State_Site), ALLOCATABLE, DIMENSION(:) :: site	! array that stores the variational parameters of the 
									! MPS variational wavefunction orgnaized in matrices 
  END TYPE Matrix_Product_State
  
  ! Structure to store the Hamiltonian, a collection of local and nearest neighbor terms
  TYPE Hamiltonian
    INTEGER lattice_length 	! lattice length
    INTEGER local_dim		! local dimension of the Hilbert space (also called 'd')
    TYPE(One_Site_Operator), ALLOCATABLE, DIMENSION(:) :: Local_Terms	! local term, i.e. operators acting on a single site
    TYPE(Two_Site_Operator), ALLOCATABLE, DIMENSION(:) :: Nearest_Neighbour_Terms	! nearest neighbor terms, i.e. operators acting on two sites
  END TYPE Hamiltonian
  
  ! Structure to store the effective Hamiltonian
  TYPE Effective_Hamiltonian
    TYPE(Matrix), ALLOCATABLE, DIMENSION(:) :: L, R	! left and right effective Hamiltonians
							! L(ii) is the effective Hamiltonian at the left of site ii (center site)
							! R(ii) is the effective Hamiltonian at the right of site ii (center site)
    TYPE(Matrix_Valued_Operator), ALLOCATABLE, DIMENSION(:) :: LC, RC	! effective Hamiltonian obtained from the contraction of a MPS on a two-site term 
									! acting on a left-center sites (LC) and center-right (RC) sites
									! the physical indices on the center site are not contracted and are the indices of the
									! of the matrix-valued operators LC and RC
									! ii is the center site of LC(ii) and RC(ii)
  END TYPE Effective_Hamiltonian
  
  ! Structure needed to describe the block structure of a two-site operator (nearest-neighbor operator)
  TYPE Map_Type
    INTEGER blk, ind 
  END TYPE Map_Type

  TYPE Tuple
    INTEGER ind1, ind2
  END TYPE Tuple
  
  TYPE Block_Type
    TYPE(Tuple), DIMENSION(:), ALLOCATABLE :: v
  END TYPE Block_Type

  TYPE Block_Info_Type
    INTEGER num_blk
    INTEGER local_dim
    INTEGER, DIMENSION(:), ALLOCATABLE :: blk_dim
    TYPE(Block_Type), DIMENSION(:), ALLOCATABLE :: block
    TYPE(Map_Type), DIMENSION(:,:), ALLOCATABLE :: map
  END TYPE Block_Info_Type

  ! Pointer to a matrix
  TYPE Matrix_Pointer 
    DOUBLE COMPLEX, DIMENSION(:,:), POINTER :: M
  END TYPE Matrix_Pointer
  
  CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Initialize_MPS_Random(d,l,m,MPS)
  ! Initialize a Matrix Product State with random numbers

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: d 	! lattice length
  INTEGER, INTENT(IN) :: l	! local dimension of the Hilbert space 
  INTEGER, INTENT(IN) :: m	! link dimension of the MPS matrices, equal for all matrices except near the edges
  TYPE(Matrix_Product_State), INTENT(OUT) :: MPS	! Matrix Product State to be initialized
  
  ! working variables
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Temp1, Temp2
  INTEGER, ALLOCATABLE, DIMENSION(:) :: LDim, RDim
  INTEGER ii, jj
  
  MPS%lattice_length = l
  MPS%local_dim = d
  !calculate consistent link dimensions
  ALLOCATE(LDim(0:MPS%lattice_length),RDim(0:MPS%lattice_length))
  LDim(0) = 1
  DO ii = 1, MPS%lattice_length
   LDim(ii) = MIN(m,LDim(ii-1)*MPS%local_dim)
  ENDDO
  RDim(l) = 1
  DO ii = MPS%lattice_length-1, 0, -1
    RDim(ii) = MIN(m,RDim(ii+1)*MPS%local_dim)
  ENDDO
  DO ii = 0, l
    LDim(ii) = MIN(LDim(ii),RDim(ii))
  ENDDO
  ! build the random MPS with the previous link dimensions
  ALLOCATE(MPS%link_dim(0:MPS%lattice_length))
  MPS%link_dim = LDim
  ALLOCATE(MPS%site(MPS%lattice_length))
  DO ii = 1, MPS%lattice_length
    IF (MPS%link_dim(ii-1)*MPS%local_dim < MPS%link_dim(ii) .OR. &
	MPS%link_dim(ii)*MPS%local_dim < MPS%link_dim(ii-1)) STOP 'Error in MPS matrix dimensions'
    ALLOCATE(MPS%site(ii)%A(MPS%local_dim))
    DO jj = 1, MPS%local_dim
      ALLOCATE(MPS%site(ii)%A(jj)%M(MPS%link_dim(ii-1),MPS%link_dim(ii)))
      ALLOCATE(Temp1(MPS%link_dim(ii-1),MPS%link_dim(ii)),Temp2(MPS%link_dim(ii-1),MPS%link_dim(ii)))
      CALL RANDOM_NUMBER(Temp1)
      CALL RANDOM_NUMBER(Temp2)
      MPS%site(ii)%A(jj)%M = Temp1 + (0.,1.)*Temp2
      DEALLOCATE(Temp1,Temp2)
    ENDDO
  ENDDO
  DEALLOCATE(LDim,RDim)
  
END SUBROUTINE Initialize_MPS_Random

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MPS_Singular_Value_Decomposition(MPS,position,Left_Right)
  ! Perform the SVD of a MPS on a site given by position and absorb the singular values on the 
  ! left (leftRight = 'L') or right (LeftRight = 'R') site of the MPS
  ! basic routine to operate on MPS, a QR decomposition would suffice instead of the SVD, not yet implemented

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: position	! site where the SVD is performed
  TYPE(Matrix_Product_State), INTENT(INOUT) :: MPS	! Input-Output: MPS to be operated on
  CHARACTER, INTENT(IN) :: Left_Right			! Flag to specify where to absorb the singular values
  
  ! working variables
  DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: M, U, V
  INTEGER dimRows, dimCols, ii
  
  IF (Left_Right == 'L' .and. position > 1 .and. position <= MPS%lattice_length) THEN
    ! SVD towards the left
    dimRows = MPS%link_dim(position-1)
    dimCols = MPS%link_dim(position)*MPS%local_dim
    ALLOCATE(M(dimRows,dimCols))
    DO ii = 1, MPS%local_dim
      M(:,1+(ii-1)*MPS%link_dim(position):ii*MPS%link_dim(position)) = MPS%site(position)%A(ii)%M
    ENDDO
    ALLOCATE(U(dimRows,dimRows))
    CALL SVD_Left(dimRows,dimCols,M,U)
    DO ii = 1, MPS%local_dim
      MPS%site(position-1)%A(ii)%M = MATMUL(MPS%site(position-1)%A(ii)%M,U)
    ENDDO
    DO ii = 1, MPS%local_dim
      MPS%site(position)%A(ii)%M = M(:,1+(ii-1)*MPS%link_dim(position):ii*MPS%link_dim(position)) 
    ENDDO
    DEALLOCATE(M,U)
  ELSE IF (Left_Right == 'R' .and. position >= 1 .and. position < MPS%lattice_length) THEN
    ! SVD towards the right
    dimRows = MPS%link_dim(position-1)*MPS%local_dim
    dimCols = MPS%link_dim(position)
    ALLOCATE(M(dimRows,dimCols))
    DO ii = 1, MPS%local_dim
      M(1+(ii-1)*MPS%link_dim(position-1):ii*MPS%link_dim(position-1),:) = MPS%site(position)%A(ii)%M
    ENDDO
    ALLOCATE(V(dimCols,dimCols))
    CALL SVD_Right(dimRows,dimCols,M,V)
    DO ii = 1, MPS%local_dim
      MPS%site(position+1)%A(ii)%M = MATMUL(V,MPS%site(position+1)%A(ii)%M)
    ENDDO
    DO ii = 1, MPS%local_dim
      MPS%site(position)%A(ii)%M = M(1+(ii-1)*MPS%link_dim(position-1):ii*MPS%link_dim(position-1),:)
    ENDDO
    DEALLOCATE(M,V)
  ENDIF
  
  CONTAINS
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  SUBROUTINE SVD_Left(dimRows,dimCols,M,U)
    ! Perform the SVD of matrix M and absorb the singular values on the left unitary matrix U
  
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dimRows,dimCols	! dimension of matrix M (dimRows x dimCols)
    DOUBLE COMPLEX, DIMENSION(dimRows,dimCols), INTENT(INOUT) :: M	! Input: matrix to be decomposed; Output: right unitary matrix
    DOUBLE COMPLEX, DIMENSION(dimRows,dimRows), INTENT(OUT)   :: U	! Output: Left unitary matrix of dimension dimRows x dimRows multiplied by the singular values

    ! working variables
    DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: V
    DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:)   :: Work
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: sigma, RWork
    INTEGER LWork, LRWork, Info, ii

    IF (dimRows > dimCols) STOP 'SVD_Left Error'
    LWork  = 24*Max(1, 2*dimRows + dimCols)
    LRWork = 5*dimRows
    ALLOCATE(Work(LWork),RWork(LRWork))
    ALLOCATE(sigma(dimRows))
    CALL ZGESVD('S','O', dimRows, dimCols, M, dimRows, sigma, U, dimRows, V, dimRows, Work, LWork, RWork, Info)
    IF (Info /= 0)  STOP 'Problems in SVD calculation'  
    DO ii = 1, dimRows
      U(:,ii) = U(:,ii)*sigma(ii)
    ENDDO
    DEALLOCATE(sigma,Work,RWork)

  END SUBROUTINE SVD_Left
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  SUBROUTINE SVD_Right(dimRows,dimCols,M,V)
    ! Perform the SVD of matrix M and absorb the singular values on the right unitary matrix V
  
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dimRows,dimCols	! dimension of matrix M (dimRows x dimCols)
    DOUBLE COMPLEX, DIMENSION(dimRows,dimCols), INTENT(INOUT) :: M	! Input: matrix to be decomposed; Output: left unitary matrix
    DOUBLE COMPLEX, DIMENSION(dimCols,dimCols), INTENT(OUT)   :: V	! Output: right unitary matrix of dimension dimCols x dimCols multiplied by the singular values
    
    ! working variables
    DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: U
    DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:)   :: Work
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: sigma, RWork
    INTEGER LWork, LRWork, Info, ii
    
    IF (dimRows < dimCols) STOP 'MPS_Singular_Value_Decomposition Error'
    LWork  = 24*Max(1, 2*dimCols + dimRows)
    LRWork = 5*dimCols
    ALLOCATE(Work(LWork),RWork(LRWork))
    ALLOCATE(sigma(dimCols))
    CALL ZGESVD('O','S', dimRows, dimCols, M, dimRows, sigma, U, dimRows, V, dimCols, Work, LWork, RWork, Info)
    IF (Info /= 0)  STOP 'Problems in SVD calculation'
    DO ii = 1, dimCols
      V(ii,:) = sigma(ii)*V(ii,:)
    ENDDO
    DEALLOCATE(sigma,Work,RWork)

  END SUBROUTINE SVD_Right
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
END SUBROUTINE MPS_Singular_Value_Decomposition

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Contract_MPS_on_Two_Site_Operator(ldim,cdim,rdim,MPS_In,op,out,LeftRight)
  ! Contract a MPS on a two-site operator, the site where the contraction occurs is the center site 
  ! If LeftRight='L' the MPS link at the left of the center site and the left indices (1) of the two-site operator are contracted
  ! If LeftRight='R' the MPS link at the right of the center site and the right indices (2) of the two-site operator are contracted
  ! The result is stored in a Matrix_Valued_Operator ('out')
  ! indexed by the physical indices on the center site. If a matrix element is not allocated
  ! the corresponding matrix elements of the two-site operator are zero

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ldim 	! dimension of the link at the left of the center site
  INTEGER, INTENT(IN) :: cdim 	! dimension of the local Hilbert space on the site to be contracted (center site)
  INTEGER, INTENT(IN) :: rdim 	! dimension of the link at the right of the center site  
  TYPE(Matrix_Product_State_Site), INTENT(IN) :: MPS_In	! MPS to be used for the contraction
  TYPE(Two_Site_Operator), INTENT(IN) :: op		! Input: two-site operator to be contracted
  TYPE(Matrix_Valued_Operator), INTENT(INOUT) :: out	! Output: matrix valued operator
  CHARACTER, INTENT(IN) :: LeftRight			! Flag to specify the left-center ('L') or center-right ('R') contraction
  
  ! working variables
  TYPE(Matrix), DIMENSION(cdim,cdim) :: Temp_Matrix1, Temp_Matrix2
  INTEGER r1, c1, r2, c2, counter, ii
  DOUBLE COMPLEX val
  
  IF (LeftRight == 'L') THEN
    counter = 0
    DO ii = 1, op%num_Elements
      r1  = op%element(ii)%r1
      c1  = op%element(ii)%c1
      r2  = op%element(ii)%r2
      c2  = op%element(ii)%c2
      val = op%element(ii)%v      
      IF (.NOT. ALLOCATED(Temp_Matrix1(r1,c1)%M)) THEN
	ALLOCATE(Temp_Matrix1(r1,c1)%M(rdim,rdim))
	Temp_Matrix1(r1,c1)%M = MATMUL(TRANSPOSE(CONJG(MPS_In%A(r1)%M)),MPS_In%A(c1)%M)
      ENDIF
      IF (.NOT. ALLOCATED(Temp_Matrix2(r2,c2)%M)) THEN
	counter = counter + 1
	ALLOCATE(Temp_Matrix2(r2,c2)%M(rdim,rdim))
	Temp_Matrix2(r2,c2)%M = val*Temp_Matrix1(r1,c1)%M
      ELSE
	Temp_Matrix2(r2,c2)%M = Temp_Matrix2(r2,c2)%M + val*Temp_Matrix1(r1,c1)%M
      ENDIF
    ENDDO
    IF (ALLOCATED(out%element)) THEN
      DO ii = 1, out%num_Elements
	DEALLOCATE(out%element(ii)%M)
      ENDDO
      DEALLOCATE(out%element)
    ENDIF
    out%num_Elements = counter
    ALLOCATE(out%element(out%num_Elements))
    counter = 0
    DO ii = 1, op%num_Elements
      r1  = op%element(ii)%r1
      c1  = op%element(ii)%c1
      r2  = op%element(ii)%r2
      c2  = op%element(ii)%c2
      IF (ALLOCATED(Temp_Matrix1(r1,c1)%M)) DEALLOCATE(Temp_Matrix1(r1,c1)%M)
      IF (ALLOCATED(Temp_Matrix2(r2,c2)%M)) THEN
	counter = counter + 1
	out%element(counter)%r = r2
	out%element(counter)%c = c2
	ALLOCATE(out%element(counter)%M(rdim,rdim))
	out%element(counter)%M = Temp_Matrix2(r2,c2)%M
	DEALLOCATE(Temp_Matrix2(r2,c2)%M)
      ENDIF
    ENDDO
    IF (counter /= out%num_Elements) STOP 'Error in Contract_MPS_on_Two_Site_Operator'
  ELSEIF (LeftRight == 'R') THEN
    counter = 0
    DO ii = 1, op%num_Elements
      r1  = op%element(ii)%r1
      c1  = op%element(ii)%c1
      r2  = op%element(ii)%r2
      c2  = op%element(ii)%c2
      val = op%element(ii)%v  
      IF (.NOT. ALLOCATED(Temp_Matrix2(r2,c2)%M)) THEN
	ALLOCATE(Temp_Matrix2(r2,c2)%M(ldim,ldim))
	Temp_Matrix2(r2,c2)%M = MATMUL(CONJG(MPS_In%A(r2)%M),TRANSPOSE(MPS_In%A(c2)%M))
      ENDIF
      IF (.NOT. ALLOCATED(Temp_Matrix1(r1,c1)%M)) THEN
	counter = counter + 1
	ALLOCATE(Temp_Matrix1(r1,c1)%M(ldim,ldim))
	Temp_Matrix1(r1,c1)%M = val*Temp_Matrix2(r2,c2)%M
      ELSE
	Temp_Matrix1(r1,c1)%M = Temp_Matrix1(r1,c1)%M + val*Temp_Matrix2(r2,c2)%M
      ENDIF
    ENDDO
    IF (ALLOCATED(out%element)) THEN
      DO ii = 1, out%num_Elements
	DEALLOCATE(out%element(ii)%M)
      ENDDO
      DEALLOCATE(out%element)
    ENDIF
    out%num_Elements = counter
    ALLOCATE(out%element(out%num_Elements))
    counter = 0
    DO ii = 1, op%num_Elements
      r1  = op%element(ii)%r1
      c1  = op%element(ii)%c1
      r2  = op%element(ii)%r2
      c2  = op%element(ii)%c2
      IF (ALLOCATED(Temp_Matrix2(r2,c2)%M)) DEALLOCATE(Temp_Matrix2(r2,c2)%M)
      IF (ALLOCATED(Temp_Matrix1(r1,c1)%M)) THEN
	counter = counter + 1
	out%element(counter)%r = r1
	out%element(counter)%c = c1
	ALLOCATE(out%element(counter)%M(ldim,ldim))
	out%element(counter)%M = Temp_Matrix1(r1,c1)%M
	DEALLOCATE(Temp_Matrix1(r1,c1)%M)
      ENDIF
    ENDDO
    IF (counter /= out%num_Elements) STOP 'Error in Contract_MPS_on_Two_Site_Operator'
  ENDIF

END SUBROUTINE Contract_MPS_on_Two_Site_Operator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Update_Effective_Hamiltonian(MPS,Ham,Eff_Ham,pos,LeftRight)
  ! Update the effective Hamiltonian using MPS and Hamiltonian 'Ham', in site 'pos'
  ! and store the result in 'Eff_Ham'
  ! Whether to update the effective Hamiltonian on the left or on the right of site 'pos'
  ! is specified by the flag LeftRight 

  IMPLICIT NONE
  INTEGER, INTENT(IN)   :: pos		! reference site 
  CHARACTER, INTENT(IN) :: LeftRight	! flags to specify whether to update the effective Hamiltonian at the left ('L') of site 'pos' 
					! or at the right ('R') of site 'pos'
  TYPE(Matrix_Product_State), INTENT(IN) :: MPS		! MPS to be used in the update
  TYPE(Hamiltonian), INTENT(IN) :: Ham			! Hamiltonian to be used in the update
  TYPE(Effective_Hamiltonian), INTENT(INOUT) :: Eff_Ham	! Output: the effective Hamiltonian
  
  ! working variables
  DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: Temp
  INTEGER r, c, counter, ii
  DOUBLE COMPLEX val
  
  IF (LeftRight == 'L' .and. pos > 1 .and. pos <= MPS%lattice_length) THEN
    ! Update at the left of site 'pos'
    ! Update  'Eff_Ham%L(pos)', namely the contraction of all Hamiltonian operators at the left (and not including) site 'pos'
    IF (ALLOCATED(Eff_Ham%L(pos)%M)) DEALLOCATE(Eff_Ham%L(pos)%M)
    ALLOCATE(Eff_Ham%L(pos)%M(MPS%link_dim(pos-1),MPS%link_dim(pos-1)),Temp(MPS%link_dim(pos-2),MPS%link_dim(pos-1)))
    Eff_Ham%L(pos)%M = 0.
    DO ii = 1, Ham%Local_Terms(pos-1)%num_Elements
      r   = Ham%Local_Terms(pos-1)%element(ii)%r
      c   = Ham%Local_Terms(pos-1)%element(ii)%c
      val = Ham%Local_Terms(pos-1)%element(ii)%v
      Eff_Ham%L(pos)%M = Eff_Ham%L(pos)%M + val*MATMUL(TRANSPOSE(CONJG(MPS%site(pos-1)%A(r)%M)),MPS%site(pos-1)%A(c)%M)
    ENDDO
    DO ii = 1, Eff_Ham%LC(pos-1)%num_Elements
      r = Eff_Ham%LC(pos-1)%element(ii)%r
      c = Eff_Ham%LC(pos-1)%element(ii)%c
      Temp = MATMUL(Eff_Ham%LC(pos-1)%element(ii)%M,MPS%site(pos-1)%A(c)%M)
      Eff_Ham%L(pos)%M = Eff_Ham%L(pos)%M + MATMUL(TRANSPOSE(CONJG(MPS%site(pos-1)%A(r)%M)),Temp)
    ENDDO
    DO ii = 1, MPS%local_dim
      Temp = MATMUL(Eff_Ham%L(pos-1)%M,MPS%site(pos-1)%A(ii)%M)
      Eff_Ham%L(pos)%M = Eff_Ham%L(pos)%M + MATMUL(TRANSPOSE(CONJG(MPS%site(pos-1)%A(ii)%M)),Temp)
    ENDDO
    DEALLOCATE(Temp)
    ! Update Eff_Ham%LC(pos), contraction of the nearest-neighbor operator acting on 'pos' and 'pos-1' 
    CALL Contract_MPS_on_Two_Site_Operator(MPS%link_dim(pos-2),MPS%local_dim,MPS%link_dim(pos-1),MPS%site(pos-1),&
					    Ham%Nearest_Neighbour_Terms(pos-1),Eff_Ham%LC(pos),'L')
  ELSEIF (LeftRight == 'R' .and. pos >= 1 .and. pos < MPS%lattice_length) THEN
    ! Update at the right of site 'pos'
    ! Update the 'Eff_Ham%R(pos)', namely the contraction of all Hamiltonian operators at the right (and not including) site 'pos'
    IF (ALLOCATED(Eff_Ham%R(pos)%M)) DEALLOCATE(Eff_Ham%R(pos)%M)
    ALLOCATE(Eff_Ham%R(pos)%M(MPS%link_dim(pos),MPS%link_dim(pos)),Temp(MPS%link_dim(pos+1),MPS%link_dim(pos)))
    Eff_Ham%R(pos)%M = 0.
    DO ii = 1, Ham%Local_Terms(pos+1)%num_Elements
      r   = Ham%Local_Terms(pos+1)%element(ii)%r
      c   = Ham%Local_Terms(pos+1)%element(ii)%c
      val = Ham%Local_Terms(pos+1)%element(ii)%v
      Eff_Ham%R(pos)%M = Eff_Ham%R(pos)%M + val*MATMUL(CONJG(MPS%site(pos+1)%A(r)%M),TRANSPOSE(MPS%site(pos+1)%A(c)%M))
    ENDDO
    DO ii = 1, Eff_Ham%RC(pos+1)%num_Elements
      r = Eff_Ham%RC(pos+1)%element(ii)%r
      c = Eff_Ham%RC(pos+1)%element(ii)%c
      Temp = MATMUL(Eff_Ham%RC(pos+1)%element(ii)%M,TRANSPOSE(MPS%site(pos+1)%A(c)%M))
      Eff_Ham%R(pos)%M = Eff_Ham%R(pos)%M + MATMUL(CONJG(MPS%site(pos+1)%A(r)%M),Temp)
    ENDDO
    DO ii = 1, MPS%local_dim
      Temp = MATMUL(Eff_Ham%R(pos+1)%M,TRANSPOSE(MPS%site(pos+1)%A(ii)%M))
      Eff_Ham%R(pos)%M = Eff_Ham%R(pos)%M + MATMUL(CONJG(MPS%site(pos+1)%A(ii)%M),Temp)
    ENDDO
    DEALLOCATE(Temp)
    ! Update Eff_Ham%RC(pos), contraction of the nearest-neighbor operator acting on 'pos' and 'pos+1' 
    CALL Contract_MPS_on_Two_Site_Operator(MPS%link_dim(pos),MPS%local_dim,MPS%link_dim(pos+1),MPS%site(pos+1),&
					    Ham%Nearest_Neighbour_Terms(pos),Eff_Ham%RC(pos),'R')
  ENDIF

END SUBROUTINE Update_Effective_Hamiltonian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MPS_Right_Gauge(MPS)
  ! Put the MPS in the right gauge, the MPS are "right unitary"

  IMPLICIT NONE
  TYPE(Matrix_Product_State), INTENT(INOUT) :: MPS
  INTEGER ii

  DO ii = MPS%lattice_length, 2, -1
    CALL MPS_Singular_Value_Decomposition(MPS,ii,'L') 
  ENDDO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MPS_Left_Gauge(MPS)
  ! Put the MPS in the left gauge, the MPS are "left unitary"

  IMPLICIT NONE
  TYPE(Matrix_Product_State), INTENT(INOUT) :: MPS
  INTEGER ii

  DO ii = 1, MPS%lattice_length-1
    CALL MPS_Singular_Value_Decomposition(MPS,ii,'R') 
  ENDDO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Normalize_MPS(MPS,LeftRight)
  ! Put a MPS in the left ('L') or right ('R') gauge and normalize

  IMPLICIT NONE
  TYPE(Matrix_Product_State), INTENT(INOUT) :: MPS
  CHARACTER, INTENT(IN) :: LeftRight
  
  ! working variables
  INTEGER ii
  DOUBLE PRECISION norm
  
  IF (LeftRight == 'L') THEN ! Left gauge
    DO ii = 1, MPS%lattice_length-1
      CALL MPS_Singular_Value_Decomposition(MPS,ii,'R')
    ENDDO
    norm  = 0.d0
    DO ii = 1, MPS%local_dim
      norm = norm + DREAL(SUM(MPS%site(MPS%lattice_length)%A(ii)%M*CONJG(MPS%site(MPS%lattice_length)%A(ii)%M)))
    ENDDO
    norm = SQRT(norm)
    DO ii = 1, MPS%local_dim
      MPS%site(MPS%lattice_length)%A(ii)%M = MPS%site(MPS%lattice_length)%A(ii)%M/norm
    ENDDO
  ELSEIF (LeftRight == 'R') THEN ! Right gauge  
    DO ii = MPS%lattice_length, 2, -1
      CALL MPS_Singular_Value_Decomposition(MPS,ii,'L')
    ENDDO
    norm  = 0.d0
    DO ii = 1, MPS%local_dim
      norm = norm + DREAL(SUM(MPS%site(1)%A(ii)%M*CONJG(MPS%site(1)%A(ii)%M)))
    ENDDO
    norm = SQRT(norm)
    DO ii = 1, MPS%local_dim
      MPS%site(1)%A(ii)%M = MPS%site(1)%A(ii)%M/norm        
    ENDDO
  ENDIF 
     
END SUBROUTINE Normalize_MPS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION Calculate_Norm_MPS(MPS,LeftRight) RESULT(norm)
  ! Put an MPS in the left ('L') or right ('R') gauge and calculate the norm

  IMPLICIT NONE
  TYPE(Matrix_Product_State), INTENT(INOUT) :: MPS
  CHARACTER, INTENT(IN) :: LeftRight
  INTEGER ii
  DOUBLE PRECISION norm
  
  IF (LeftRight == 'L') THEN ! Left gauge
    DO ii = 1, MPS%lattice_length-1
      CALL MPS_Singular_Value_Decomposition(MPS,ii,'R')
    ENDDO
    norm  = 0.d0
    DO ii = 1, MPS%local_dim
      norm = norm + DREAL(SUM(MPS%site(MPS%lattice_length)%A(ii)%M*CONJG(MPS%site(MPS%lattice_length)%A(ii)%M)))
    ENDDO
    norm = SQRT(norm)
  ELSEIF (LeftRight == 'R') THEN ! Right gauge
    DO ii = MPS%lattice_length, 2, -1
      CALL MPS_Singular_Value_Decomposition(MPS,ii,'L')
    ENDDO
    norm  = 0.d0
    DO ii = 1, MPS%local_dim
      norm = norm + DREAL(SUM(MPS%site(1)%A(ii)%M*CONJG(MPS%site(1)%A(ii)%M)))
    ENDDO
    norm = SQRT(norm)
  ENDIF
          
END FUNCTION Calculate_Norm_MPS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
SUBROUTINE Initialize_Effective_Hamiltonian(MPS,Ham,Eff_Ham)
  ! Initialize the effective Hamiltonian

  IMPLICIT NONE
  TYPE(Matrix_Product_State), INTENT(INOUT) :: MPS
  TYPE(Hamiltonian), INTENT(IN) :: Ham
  TYPE(Effective_Hamiltonian), INTENT(OUT) :: Eff_Ham
  INTEGER ii, jj

  CALL MPS_Right_Gauge(MPS)
  ALLOCATE(Eff_Ham%L(MPS%lattice_length),Eff_Ham%R(MPS%lattice_length),&
	   Eff_Ham%LC(MPS%lattice_length),Eff_Ham%RC(MPS%lattice_length))
  ALLOCATE(Eff_Ham%L(1)%M(1,1))
  Eff_Ham%L(1)%M(1,1) = 0.
  EFF_Ham%LC(1)%num_Elements = 0
  ALLOCATE(Eff_Ham%R(Ham%lattice_length)%M(1,1))
  Eff_Ham%R(Ham%lattice_length)%M(1,1) = 0.
  Eff_Ham%RC(Ham%lattice_length)%num_Elements = 0
  DO ii = Ham%lattice_length-1, 1, -1
    CALL Update_Effective_Hamiltonian(MPS,Ham,Eff_Ham,ii,'R')
  ENDDO

END SUBROUTINE Initialize_Effective_Hamiltonian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SUBROUTINE Print_Effective_Hamiltonian_R(Eff_Ham,pos)
!   
!   IMPLICIT NONE
!   INTEGER, INTENT(IN) :: pos
!   TYPE(Effective_Hamiltonian), INTENT(IN) :: Eff_Ham
!   INTEGER ii, jj
!   INTEGER, DIMENSION(2) :: m
!   
!   20 FORMAT(20('   ',ES10.3,',',ES10.3,''))
!   print *, '------------------------------------------------------------------------------------------'
!   m = SHAPE(Eff_Ham%R(pos)%M)
!   DO ii = 1, m(1)
!     print 20, Eff_Ham%R(pos)%M(ii,:)
!   ENDDO
!   print *, '------------------------------------------------------------------------------------------'
!   DO ii = 1, Eff_Ham%RC(pos)%num_Elements
!     print *, Eff_Ham%RC(pos)%element(ii)%r, Eff_Ham%RC(pos)%element(ii)%c
!     m = SHAPE(Eff_Ham%RC(pos)%element(ii)%M)
!     DO jj = 1, m(1)
!       print 20, Eff_Ham%RC(pos)%element(ii)%M(jj,:)
!     ENDDO
!   ENDDO
!   print *, '------------------------------------------------------------------------------------------'
! 
!   
! END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SUBROUTINE Print_Effective_Hamiltonian_L(Eff_Ham,pos)
!   
!   IMPLICIT NONE
!   INTEGER, INTENT(IN) :: pos
!   TYPE(Effective_Hamiltonian), INTENT(IN) :: Eff_Ham  
!   INTEGER ii, jj
!   INTEGER, DIMENSION(2) :: m
!   
!   20 FORMAT(20('   ',ES10.3,',',ES10.3,''))
!   print *, '------------------------------------------------------------------------------------------'
!   m = SHAPE(Eff_Ham%L(pos)%M)
!   DO ii = 1, m(1)
!     print 20, Eff_Ham%L(pos)%M(ii,:)
!   ENDDO
!   print *, '------------------------------------------------------------------------------------------'
!   DO ii = 1, Eff_Ham%LC(pos)%num_Elements
!     print *, Eff_Ham%LC(pos)%element(ii)%r, Eff_Ham%LC(pos)%element(ii)%c
!     m = SHAPE(Eff_Ham%LC(pos)%element(ii)%M)
!     DO jj = 1, m(1)
!       print 20, Eff_Ham%LC(pos)%element(ii)%M(jj,:)
!     ENDDO
!   ENDDO
!   print *, '------------------------------------------------------------------------------------------'
! 
!   
! END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Multiply_Ham_Psi_One_Site(pos,ldim,rdim,cdim,Psi_In,Psi_Out,Ham,Eff_Ham)
  ! Multiply the effective Hamiltonian 'Eff_Ham' relative to a single site by a vector 'Psi_In'
  ! and store the result in 'Psi_out'
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: pos	! site where the effective Hamiltonian has been calculated
  INTEGER, INTENT(IN) :: ldim 	! dimension of the link at the left of the center site 'pos'
  INTEGER, INTENT(IN) :: cdim 	! dimension of the local Hilbert space on the center site 'pos'
  INTEGER, INTENT(IN) :: rdim 	! dimension of the link at the right of the center site 'pos'
  DOUBLE COMPLEX, DIMENSION(ldim,rdim,cdim), INTENT(IN)  :: Psi_In	! input vector
  DOUBLE COMPLEX, DIMENSION(ldim,rdim,cdim), INTENT(OUT) :: Psi_Out	! output vector
  TYPE(Hamiltonian), INTENT(IN) :: Ham					! Hamiltonian
  TYPE(Effective_Hamiltonian), INTENT(IN) :: Eff_Ham			! effective Hamiltonian (contracted Hamiltonian)

  ! working variables
  INTEGER r, c, ii
  DOUBLE COMPLEX val
  
  Psi_Out = 0.
  DO ii = 1, Ham%Local_Terms(pos)%num_Elements  ! Local operator
    r   = Ham%Local_Terms(pos)%element(ii)%r
    c   = Ham%Local_Terms(pos)%element(ii)%c
    val = Ham%Local_Terms(pos)%element(ii)%v
    Psi_Out(:,:,r) = Psi_Out(:,:,r) + val*Psi_In(:,:,c)
  ENDDO
  DO ii = 1, cdim  			! Effective Hamiltonian Left and Right Block 
    Psi_Out(:,:,ii) = Psi_Out(:,:,ii) + MATMUL(Eff_Ham%L(pos)%M,Psi_In(:,:,ii)) &
				      + MATMUL(Psi_In(:,:,ii),TRANSPOSE(Eff_Ham%R(pos)%M))
  ENDDO
  DO ii = 1, Eff_Ham%LC(pos)%num_Elements	! Effective Hamiltonian Left-Center Block 
      r = Eff_Ham%LC(pos)%element(ii)%r
      c = Eff_Ham%LC(pos)%element(ii)%c
      Psi_Out(:,:,r) = Psi_Out(:,:,r) + MATMUL(Eff_Ham%LC(pos)%element(ii)%M,Psi_In(:,:,c))
  ENDDO
  DO ii = 1, Eff_Ham%RC(pos)%num_Elements	! Effective Hamiltonian Right-Center Block
      r = Eff_Ham%RC(pos)%element(ii)%r
      c = Eff_Ham%RC(pos)%element(ii)%c
      Psi_Out(:,:,r) = Psi_Out(:,:,r) + MATMUL(Psi_In(:,:,c),TRANSPOSE(Eff_Ham%RC(pos)%element(ii)%M))
  ENDDO

END SUBROUTINE Multiply_Ham_Psi_One_Site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Multiply_Ham_Psi_Two_Site(pos,ldim,rdim,cdim,Psi_In,Psi_Out,Ham,Eff_Ham)
  ! Multiply the effective Hamiltonian 'Eff_Ham' relative to a two sites by a vector 'Psi_In'
  ! and store the result in 'Psi_out'
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: pos	! site where the effective Hamiltonian has been calculated (two sites 'pos' and 'pos+1')
  INTEGER, INTENT(IN) :: ldim 	! dimension of the link at the left of the center site 'pos'
  INTEGER, INTENT(IN) :: cdim 	! dimension of the local Hilbert space on the center site 'pos'
  INTEGER, INTENT(IN) :: rdim 	! dimension of the link at the right of the site 'pos+1'
  DOUBLE COMPLEX, DIMENSION(ldim,rdim,cdim,cdim), INTENT(IN)  :: Psi_In		! input vector
  DOUBLE COMPLEX, DIMENSION(ldim,rdim,cdim,cdim), INTENT(OUT) :: Psi_Out	! output vector
  TYPE(Hamiltonian), INTENT(IN) :: Ham						! Hamiltonian
  TYPE(Effective_Hamiltonian), INTENT(IN) :: Eff_Ham				! effective Hamiltonian (contracted Hamiltonian)
  
  ! working variables
  INTEGER r1, c1, r2, c2, ii, jj
  DOUBLE COMPLEX val
  
  Psi_Out = 0.
  DO ii = 1, Ham%Local_Terms(pos)%num_Elements		! Local operator on site 'pos'
    r1  = Ham%Local_Terms(pos)%element(ii)%r
    c1  = Ham%Local_Terms(pos)%element(ii)%c
    val = Ham%Local_Terms(pos)%element(ii)%v
    DO jj = 1, cdim
      Psi_Out(:,:,r1,jj) = Psi_Out(:,:,r1,jj) + val*Psi_In(:,:,c1,jj)
    ENDDO
  ENDDO
  DO ii = 1, Ham%Local_Terms(pos+1)%num_Elements 	! Local operator on site 'pos+1'
    r2  = Ham%Local_Terms(pos+1)%element(ii)%r
    c2  = Ham%Local_Terms(pos+1)%element(ii)%c
    val = Ham%Local_Terms(pos+1)%element(ii)%v
    DO jj = 1, cdim
      Psi_Out(:,:,jj,r2) = Psi_Out(:,:,jj,r2) + val*Psi_In(:,:,jj,c2)
    ENDDO
  ENDDO
  DO ii = 1, Ham%Nearest_Neighbour_Terms(pos)%num_Elements	! Nearest-neighbor terms on sites 'pos' and 'pos+1'
    r1  = Ham%Nearest_Neighbour_Terms(pos)%element(ii)%r1
    c1  = Ham%Nearest_Neighbour_Terms(pos)%element(ii)%c1
    r2  = Ham%Nearest_Neighbour_Terms(pos)%element(ii)%r2
    c2  = Ham%Nearest_Neighbour_Terms(pos)%element(ii)%c2
    val = Ham%Nearest_Neighbour_Terms(pos)%element(ii)%v
    Psi_Out(:,:,r1,r2) = Psi_Out(:,:,r1,r2) + val*Psi_In(:,:,c1,c2)
  ENDDO
  DO ii = 1, cdim 			! Effective Hamiltonian Left and Right Block 
    DO jj = 1, cdim
      Psi_Out(:,:,ii,jj) = Psi_Out(:,:,ii,jj) + MATMUL(Eff_Ham%L(pos)%M,Psi_In(:,:,ii,jj)) &
					      + MATMUL(Psi_In(:,:,ii,jj),TRANSPOSE(Eff_Ham%R(pos+1)%M))
    ENDDO
  ENDDO
  DO ii = 1, Eff_Ham%LC(pos)%num_Elements	! Effective Hamiltonian Left-Center Block 
    r1 = Eff_Ham%LC(pos)%element(ii)%r
    c1 = Eff_Ham%LC(pos)%element(ii)%c
    DO jj = 1, cdim
      Psi_Out(:,:,r1,jj) = Psi_Out(:,:,r1,jj) + MATMUL(Eff_Ham%LC(pos)%element(ii)%M,Psi_In(:,:,c1,jj))
    ENDDO
  ENDDO
  DO ii = 1, Eff_Ham%RC(pos+1)%num_Elements	! Effective Hamiltonian Center-Right Block
    r2 = Eff_Ham%RC(pos+1)%element(ii)%r
    c2 = Eff_Ham%RC(pos+1)%element(ii)%c
    DO jj = 1, cdim
      Psi_Out(:,:,jj,r2) = Psi_Out(:,:,jj,r2) + MATMUL(Psi_In(:,:,jj,c2),TRANSPOSE(Eff_Ham%RC(pos+1)%element(ii)%M))
    ENDDO
  ENDDO

END SUBROUTINE Multiply_Ham_Psi_Two_Site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Build_Hamiltonian_Matrix(n,H,pos,ldim,rdim,cdim,Ham,Eff_Ham)
  ! build the effective Hamiltonian in full matrix form
  ! used to test the Davidson algorithm for small MPS link dimension (Test routine)
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, pos, ldim, rdim, cdim
  DOUBLE COMPLEX, DIMENSION(n,n), INTENT(OUT) :: H
  TYPE(Hamiltonian), INTENT(IN) :: Ham
  TYPE(Effective_Hamiltonian), INTENT(IN) :: Eff_Ham
  INTEGER ii, jj, kk, mm, nn, rowI, colI, el, r, c
  DOUBLE COMPLEX val
  
  H = 0.
  DO el = 1, Ham%Local_Terms(pos)%num_Elements
    r   = Ham%Local_Terms(pos)%element(el)%r
    c   = Ham%Local_Terms(pos)%element(el)%c
    val = Ham%Local_Terms(pos)%element(el)%v
    DO mm = 1, rdim
      DO ii = 1, ldim
	rowI = Map(ii,ldim,mm,rdim,r)
	colI = Map(ii,ldim,mm,rdim,c)
	H(rowI,colI) = H(rowI,colI) + val 
      ENDDO
    ENDDO
  ENDDO
  DO kk = 1, cdim
    DO mm = 1, rdim
      DO ii = 1, ldim
	DO jj = 1, ldim
	  rowI = Map(ii,ldim,mm,rdim,kk)
	  colI = Map(jj,ldim,mm,rdim,kk)
	  H(rowI,colI) = H(rowI,colI) + Eff_Ham%L(pos)%M(ii,jj)
	ENDDO
      ENDDO
    ENDDO
  ENDDO
  DO kk = 1, cdim
    DO mm = 1, rdim
      DO nn = 1, rdim
	DO ii = 1, ldim
	  rowI = Map(ii,ldim,mm,rdim,kk)
	  colI = Map(ii,ldim,nn,rdim,kk)
	  H(rowI,colI) = H(rowI,colI) + Eff_Ham%R(pos)%M(mm,nn)
	ENDDO
      ENDDO
    ENDDO
  ENDDO
  DO el = 1, Eff_Ham%LC(pos)%num_Elements
    r = Eff_Ham%LC(pos)%element(el)%r
    c = Eff_Ham%LC(pos)%element(el)%c
    DO mm = 1, rdim
      DO ii = 1, ldim
	DO jj = 1, ldim
	  rowI = Map(ii,ldim,mm,rdim,r)
	  colI = Map(jj,ldim,mm,rdim,c)
	  H(rowI,colI) = H(rowI,colI) + Eff_Ham%LC(pos)%element(el)%M(ii,jj)
	ENDDO
      ENDDO
    ENDDO
  ENDDO
  DO el = 1, Eff_Ham%RC(pos)%num_Elements
    r = Eff_Ham%RC(pos)%element(el)%r
    c = Eff_Ham%RC(pos)%element(el)%c
    DO mm = 1, rdim
      DO nn = 1, rdim
	DO ii = 1, ldim
	  rowI = Map(ii,ldim,mm,rdim,r)
	  colI = Map(ii,ldim,nn,rdim,c)
	  H(rowI,colI) = H(rowI,colI) + Eff_Ham%RC(pos)%element(el)%M(mm,nn)
	ENDDO
      ENDDO
    ENDDO
  ENDDO
  
  CONTAINS
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  PURE FUNCTION Map(ii,ldim,mm,rdim,kk)
    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ii, ldim, mm, rdim, kk
    INTEGER Map
    Map = ii + ldim*(mm-1) + ldim*rdim*(kk-1)
    
  END FUNCTION Map
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
END SUBROUTINE Build_Hamiltonian_Matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Tensor3_to_Array(dim1, dim2, dim3, tensor, array_dim, array)
  ! Map a three-index array (Tensor3) to an array

  INTEGER, INTENT(IN) :: dim1, dim2, dim3, array_dim
  DOUBLE COMPLEX, DIMENSION(dim1, dim2, dim3), INTENT(IN) :: tensor
  DOUBLE COMPLEX, DIMENSION(array_dim), INTENT(OUT) :: array
  INTEGER ii, jj, kk, counter
  
  counter = 0
  DO kk = 1, dim3
    DO jj = 1, dim2
      DO ii = 1, dim1
	counter = counter + 1
	array(counter) = tensor(ii,jj,kk)
      ENDDO
    ENDDO
  ENDDO
  
END SUBROUTINE Tensor3_to_Array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Array_to_Tensor3(dim1, dim2, dim3, tensor, array_dim, array)
  ! Map an array to a three-index array (Tensor3)

  INTEGER, INTENT(IN) :: dim1, dim2, dim3, array_dim
  DOUBLE COMPLEX, DIMENSION(dim1, dim2, dim3), INTENT(OUT) :: tensor
  DOUBLE COMPLEX, DIMENSION(array_dim), INTENT(IN) :: array
  INTEGER ii, jj, kk, counter
  
  counter = 0
  DO kk = 1, dim3
    DO jj = 1, dim2
      DO ii = 1, dim1
	counter = counter + 1
	tensor(ii,jj,kk) = array(counter)
      ENDDO
    ENDDO
  ENDDO
  
END SUBROUTINE Array_to_Tensor3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Tensor4_to_Array(dim1, dim2, dim3, dim4, tensor, array_dim, array)
  ! Map a four-index array (Tensor4) to an array

  INTEGER, INTENT(IN) :: dim1, dim2, dim3, dim4, array_dim
  DOUBLE COMPLEX, DIMENSION(dim1, dim2, dim3,dim4), INTENT(IN) :: tensor
  DOUBLE COMPLEX, DIMENSION(array_dim), INTENT(OUT) :: array
  INTEGER ii, jj, kk, ll, counter
  
  counter = 0
  DO ll = 1, dim4
    DO kk = 1, dim3
      DO jj = 1, dim2
	DO ii = 1, dim1
	  counter = counter + 1
	  array(counter) = tensor(ii,jj,kk,ll)
	ENDDO
      ENDDO
    ENDDO
  ENDDO
  
END SUBROUTINE Tensor4_to_Array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Array_to_Tensor4(dim1, dim2, dim3, dim4, tensor, array_dim, array)
  ! Map an array to a four-index array (Tensor4)

  INTEGER, INTENT(IN) :: dim1, dim2, dim3, dim4, array_dim
  DOUBLE COMPLEX, DIMENSION(dim1, dim2, dim3,dim4), INTENT(OUT) :: tensor
  DOUBLE COMPLEX, DIMENSION(array_dim), INTENT(IN) :: array
  INTEGER ii, jj, kk, ll, counter
  
  counter = 0
  DO ll = 1, dim4
    DO kk = 1, dim3
      DO jj = 1, dim2
	DO ii = 1, dim1
	  counter = counter + 1
	  tensor(ii,jj,kk,ll) = array(counter)
	ENDDO
      ENDDO
    ENDDO
  ENDDO
  
END SUBROUTINE Array_to_Tensor4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Davidson(eigenvalue,n,eigenvector)
  ! wrapper subroutine for the davidson diagonalization

  IMPLICIT NONE
  INTEGER,     	    INTENT(IN)  :: n		! dimension of 'eigenvector'
  DOUBLE PRECISION, INTENT(OUT) :: eigenvalue	! the eigenvalue obtained form the Davidson diagonalization 
  DOUBLE COMPLEX,  DIMENSION(n) :: eigenvector	! the eigenvector (ground state wavefunction)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  LOGICAL		  :: use_Guess = .TRUE. ! If use_Guess=.true. use the current MPS matrices to start iteration
  LOGICAL 	          :: wanted    = .TRUE. ! If wanted=.true. then computes the converged eigenvectors
  INTEGER                 :: kmax,jmax,jmin,maxstep,method,m,l,maxnmv,order,testspace,lwork
  INTEGER 		  :: KmaxUser = 2
  DOUBLE PRECISION        :: tol,lock,targetEn,norm,emin,etemp
  DOUBLE COMPLEX,   DIMENSION(:),   ALLOCATABLE :: alpha,beta,tmp,residu
  DOUBLE COMPLEX,   DIMENSION(:,:), ALLOCATABLE :: eivec,zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  kmax = 1
  targetEn = 0.
  tol = 1.d-12    ! Tolerance of the eigensolutions: $\Vert \beta H_{SB} x - \alpha x \vert$
  maxnmv = 100    ! Maximum number of matvecs in cgstab or gmres (very problem dependent; typically choose in [5-100])
  order = -1      ! Selection criterion for Ritz values:  0 (nearest to target);  -1 (smallest real part)
  IF(order == 0)  testspace = 3 ! put 3 if a reasonable value for target is known, else take 2
  IF(order /= 0)  testspace = 2

  IF (3*KmaxUser <= 20) jmax=20          ! Maximum size of the search space:
  IF (3*KmaxUser >  20) jmax=3*KmaxUser
  jmin=2*KmaxUser                        ! Minimum size of the search space
  maxstep = 1000                         ! Maximum number of Jacobi-Davidson iterations
  lock = 1.0d-12                         ! Tracking parameter
  method = 2                             ! Method for the linear equation solver  1: gmres(m)  2: cgstab(l)
  m = 30                                 ! Maximum dimension of searchspace for gmres(m):
  l= 2                                   ! Degree of gmres-polynomial in bi-cgstab(l):
  IF (method == 1) lwork =  4 +  m  + 5*jmax + 3*KmaxUser  !Size of workspace
  IF (method == 2) lwork = 10 + 6*l + 5*jmax + 3*KmaxUser  !KmaxUser is used since Kmax = 1 gives problems ...!

  ALLOCATE (alpha(jmax),beta(jmax),eivec(n,kmax))
  alpha =0.d0
  beta  =0.d0
  eiVec =0.d0
  ALLOCATE (tmp(n),residu(n),zwork(n,lwork))
  tmp   =0.d0
  residu=0.d0
  zwork =0.d0
  
  CALL JDQZ(alpha, beta, EIVEC, wanted, n, targetEn, tol, Kmax, jmax, jmin, method, m, l, maxnmv, maxstep, &
            lock, order, testspace, zwork, lwork, use_Guess)

  CALL AMUL  ( n, eivec(1,1), residu )
  CALL ZSCAL ( n, beta(1), residu, 1 )
  CALL BMUL  ( n, eivec(1,1), tmp )
  CALL ZAXPY ( n, -alpha(1), tmp, 1, residu, 1 )
  DEALLOCATE (zwork,tmp,residu)

  eigenvalue = DBLE(alpha(1)/beta(1))
  DEALLOCATE (alpha,beta)

  eigenvector=eivec(:,1)
  norm=SUM(CONJG(eigenvector)*eigenvector)
  eigenvector=eigenvector/(norm**0.5d0)
  DEALLOCATE (eivec)

END SUBROUTINE Davidson

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MPS_Multiply_and_Truncate(MPS,op,pos,m_Min,m_Max,eps,LeftRight,coeff)
  ! Multiply a MPS by a two-site operator and truncate after the SVD
  ! Used in White correction and for the time evolution

  IMPLICIT NONE
  CHARACTER, INTENT(IN) :: LeftRight		! Flag to specify where to absorb the singular values
  INTEGER, INTENT(IN) 	:: pos			! multiply the MPS on sites 'pos' and 'pos+1' 
  INTEGER, INTENT(IN) :: m_Min			! minimum link dimension
  INTEGER, INTENT(IN) :: m_Max			! maximum link dimension 
  DOUBLE PRECISION, INTENT(IN) :: eps		! truncation error
  DOUBLE PRECISION, INTENT(IN) :: coeff		! multiplication coefficient
  TYPE(Matrix_Product_State), INTENT(INOUT) :: MPS	! Matrix Product State
  TYPE(Two_Site_Operator), INTENT(IN) :: op		! multiplication operator
  
  ! working variables
  DOUBLE COMPLEX, DIMENSION(:,:), POINTER :: B
  TYPE(Matrix), DIMENSION(MPS%local_dim,MPS%local_dim) :: Temp
  DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE :: U, V
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: S
  INTEGER ii, jj, d, ld, rd, el, dimMin, r1, c1, r2, c2
  DOUBLE COMPLEX val
  DOUBLE PRECISION ErrorSum
  
  TYPE Matrix_Pointer 
    DOUBLE COMPLEX, DIMENSION(:,:), POINTER :: M
  END TYPE Matrix_Pointer
  TYPE(Matrix_Pointer), DIMENSION(MPS%local_dim,MPS%local_dim) :: BP  
  
  ld = MPS%link_dim(pos-1)
  rd = MPS%link_dim(pos+1)
  d = MPS%local_dim
  ! Multiplication and contraction of middle link
  ALLOCATE(B(ld*d,rd*d))
  DO ii = 1, d
    DO jj = 1, d 
      BP(ii,jj)%M => B(1+(ii-1)*ld:ii*ld,1+(jj-1)*rd:jj*rd) 
    ENDDO
  ENDDO 
  B = 0.
  DO el = 1, op%num_Elements
    r1  = op%element(el)%r1
    c1  = op%element(el)%c1
    r2  = op%element(el)%r2
    c2  = op%element(el)%c2
    val = op%element(el)%v
    IF (.NOT. ALLOCATED(Temp(c1,c2)%M)) THEN
      ALLOCATE(Temp(c1,c2)%M(ld,rd))
      Temp(c1,c2)%M = MATMUL(MPS%site(pos)%A(c1)%M,MPS%site(pos+1)%A(c2)%M)
    ENDIF
    BP(r1,r2)%M = BP(r1,r2)%M + coeff*val*Temp(c1,c2)%M 
  ENDDO
  ! SVD
  dimMin = MIN(ld,rd)*d
  ALLOCATE(U(ld*d,dimMin),V(dimMin,rd*d),S(dimMin))
  CALL SVD(ld*d,dimMin,rd*d,B,U,S,V)
  ! Truncation
  S = S/SQRT(SUM(S**2))
  ErrorSum = 0.
  IF (dimMin > m_Min) THEN
    loop1: DO ii = dimMin, 1, -1
        IF ( ii <= m_Max .AND. ((ErrorSum + S(ii)**2) > eps) .OR. ii <= m_Min) THEN
           MPS%link_dim(pos) = ii
           EXIT loop1
        ENDIF
        ErrorSum = ErrorSum + S(ii)**2
    ENDDO loop1
  ENDIF
  IF (LeftRight == 'L') THEN
    DO ii = 1, MPS%link_dim(pos)
      U(:,ii) = U(:,ii)*S(ii) 
    ENDDO
  ELSEIF (LeftRight == 'R') THEN
    DO ii = 1, MPS%link_dim(pos)
      V(ii,:) = S(ii)*V(ii,:) 
    ENDDO  
  ENDIF
  DO ii = 1, d
    DEALLOCATE(MPS%site(pos)%A(ii)%M)
    DEALLOCATE(MPS%site(pos+1)%A(ii)%M)
    ALLOCATE(MPS%site(pos)%A(ii)%M(ld,MPS%link_dim(pos)))
    ALLOCATE(MPS%site(pos+1)%A(ii)%M(MPS%link_dim(pos),rd))
    MPS%site(pos)%A(ii)%M = U(1+(ii-1)*ld:ii*ld,1:MPS%link_dim(pos))
    MPS%site(pos+1)%A(ii)%M = V(1:MPS%link_dim(pos),1+(ii-1)*rd:ii*rd)
  ENDDO
  DO el = 1, op%num_Elements
    c1 = op%element(el)%c1
    c2 = op%element(el)%c2
    IF (ALLOCATED(Temp(c1,c2)%M)) DEALLOCATE(Temp(c1,c2)%M)
  ENDDO
  DEALLOCATE(B,U,S,V)

  CONTAINS
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  SUBROUTINE SVD(dimRow,dimMin,dimCol,M,U,S,V)
  ! SVD of a matrix M of dimension dimRows x dimCols
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: dimRow, dimCol, dimMin
  DOUBLE COMPLEX, DIMENSION(dimRow,dimCol), INTENT(IN)  :: M
  DOUBLE COMPLEX, DIMENSION(dimRow,dimMin), INTENT(OUT) :: U
  DOUBLE COMPLEX, DIMENSION(dimMin,dimCol), INTENT(OUT) :: V
  DOUBLE PRECISION, DIMENSION(dimMin), INTENT(OUT) :: S
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
  
END SUBROUTINE MPS_Multiply_and_Truncate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MPS_Multiply_and_Truncate_with_Particle_Number_Correction(MPS,op,pos,m_Min,m_Max,eps,LeftRight,correction)
  ! Multiply a MPS by a two-site operator and truncate after the SVD
  ! used for the time evolution. The modified truncation procedure can be used in order to conserve the particle number


  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: correction		! flag to specify if the correction is used or not
  CHARACTER, INTENT(IN) :: LeftRight		! Flag to specify where to absorb the singular values
  INTEGER, INTENT(IN) 	:: pos			! multiply the MPS on sites 'pos' and 'pos+1' 
  INTEGER, INTENT(IN) :: m_Min			! minimum link dimension
  INTEGER, INTENT(IN) :: m_Max			! maximum link dimension 
  DOUBLE PRECISION, INTENT(IN) :: eps		! truncation error
  TYPE(Matrix_Product_State), INTENT(INOUT) :: MPS	! Matrix Product State
  TYPE(Two_Site_Operator), INTENT(IN) :: op		! multiplication operator

  ! working variables
  DOUBLE COMPLEX, DIMENSION(:,:), POINTER :: B, C
  TYPE(Matrix), DIMENSION(MPS%local_dim,MPS%local_dim) :: Temp
  DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE :: U, V
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: S
  INTEGER ii, jj, d, ld, rd, el, dimMin, r1, c1, r2, c2
  DOUBLE COMPLEX val
  DOUBLE PRECISION ErrorSum
  
  TYPE(Matrix_Pointer), DIMENSION(MPS%local_dim,MPS%local_dim) :: BP  
  
  ld = MPS%link_dim(pos-1)
  rd = MPS%link_dim(pos+1)
  d = MPS%local_dim
  ! Multiplication and contraction of middle link
  ALLOCATE(B(ld*d,rd*d),C(ld*d,rd*d))
  DO ii = 1, d
    DO jj = 1, d 
      BP(ii,jj)%M => B(1+(ii-1)*ld:ii*ld,1+(jj-1)*rd:jj*rd) 
    ENDDO
  ENDDO
  B = 0.
  DO el = 1, op%num_Elements
    r1  = op%element(el)%r1
    c1  = op%element(el)%c1
    r2  = op%element(el)%r2
    c2  = op%element(el)%c2
    val = op%element(el)%v
    IF (.NOT. ALLOCATED(Temp(c1,c2)%M)) THEN
      ALLOCATE(Temp(c1,c2)%M(ld,rd))
      Temp(c1,c2)%M = MATMUL(MPS%site(pos)%A(c1)%M,MPS%site(pos+1)%A(c2)%M)
    ENDIF
    BP(r1,r2)%M = BP(r1,r2)%M + val*Temp(c1,c2)%M 
  ENDDO
  DO el = 1, op%num_Elements
    c1 = op%element(el)%c1
    c2 = op%element(el)%c2
    IF (ALLOCATED(Temp(c1,c2)%M)) DEALLOCATE(Temp(c1,c2)%M)
  ENDDO  
  dimMin = MIN(ld,rd)*d
  ALLOCATE(U(ld*d,dimMin),V(dimMin,rd*d),S(dimMin))
  C = B
  CALL SVD(ld*d,dimMin,rd*d,C,U,S,V)
  ! Truncation
  S = S/SQRT(SUM(S**2))
  ErrorSum = 0.d0
  IF (dimMin > m_Min) THEN
    loop1: DO ii = dimMin, 1, -1
        IF ( ii <= m_Max .AND. ((ErrorSum + S(ii)**2) > eps) .OR. ii <= m_Min) THEN
           MPS%link_dim(pos) = ii
           EXIT loop1
        ENDIF
        ErrorSum = ErrorSum + S(ii)**2
    ENDDO loop1
  ENDIF
  
  ! Absorb the singular values equally on the left and right unitary matrices
  ! Suitable initial condition for particle number correction algorithm
  S(:MPS%link_dim(pos)) = S(:MPS%link_dim(pos))/SQRT(SUM(S(:MPS%link_dim(pos))**2))
  DO jj = 1, MPS%link_dim(pos)
     U(:,jj) = U(:,jj)*SQRT(S(jj)) 
     V(jj,:) = SQRT(S(jj))*V(jj,:) 
  ENDDO
  ! Store in the MPS
  DO ii = 1, d
    DEALLOCATE(MPS%site(pos)%A(ii)%M)
    DEALLOCATE(MPS%site(pos+1)%A(ii)%M)
    ALLOCATE(MPS%site(pos)%A(ii)%M(ld,MPS%link_dim(pos)))
    ALLOCATE(MPS%site(pos+1)%A(ii)%M(MPS%link_dim(pos),rd))
    MPS%site(pos)%A(ii)%M = U(1+(ii-1)*ld:ii*ld,1:MPS%link_dim(pos))
    MPS%site(pos+1)%A(ii)%M = V(1:MPS%link_dim(pos),1+(ii-1)*rd:ii*rd)
  ENDDO
  DEALLOCATE(U,V,C)

  IF (correction) CALL Particle_Number_Conservation_Correction(MPS,pos,S,BP)
!  print *, 'Corrected'
  
  ! After the particle number correction perform one more SVD to restore the gauge 
  ! on the corrected matrices
  IF (LeftRight == 'L') THEN
    CALL MPS_Singular_Value_Decomposition(MPS,pos+1,LeftRight)  
  ELSEIF (LeftRight == 'R') THEN
    CALL MPS_Singular_Value_Decomposition(MPS,pos,LeftRight)
  ENDIF
  
  DEALLOCATE(B,S)
  
  CONTAINS
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  SUBROUTINE SVD(dimRow,dimMin,dimCol,M,U,S,V)
  ! SVD of a matrix M of dimension dimRows x dimCols
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: dimRow, dimCol, dimMin			
  DOUBLE COMPLEX, DIMENSION(dimRow,dimCol), INTENT(IN)  :: M
  DOUBLE COMPLEX, DIMENSION(dimRow,dimMin), INTENT(OUT) :: U
  DOUBLE COMPLEX, DIMENSION(dimMin,dimCol), INTENT(OUT) :: V
  DOUBLE PRECISION, DIMENSION(dimMin), INTENT(OUT) :: S
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
  
END SUBROUTINE MPS_Multiply_and_Truncate_with_Particle_Number_Correction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Particle_Number_Conservation_Correction(MPS,pos,S,BP)
  ! This is the routine where the algorithm described in the paper is implemented

  IMPLICIT NONE
  TYPE(Matrix_Product_State), INTENT(INOUT) :: MPS	! Input: The MPS with the result of the SVD and truncation already performed
							! singular values are absorbed evenly between matrices in site 'pos' and 'pos+1'
							! Output: The MPS with corrected truncated matrices
  INTEGER, INTENT(IN) :: pos				! the procedures corrects the MPS matrices on site 'pos' and 'pos+1'
  DOUBLE PRECISION, DIMENSION(MPS%link_dim(pos)), INTENT(IN) :: S	! singular values from the SVD and already truncated 
  TYPE(Matrix_Pointer), DIMENSION(MPS%local_dim,MPS%local_dim), INTENT(IN) :: BP	! Full time-evolved MPS before SVD and truncation
  
  ! working variables
  TYPE(Matrix), DIMENSION(MPS%local_dim) :: A1, A2, A1_prev, A2_prev, M1, M2
  DOUBLE COMPLEX, DIMENSION(MPS%link_dim(pos),MPS%link_dim(pos)) :: C_id, C_n, C_nn, D_id, D_n, D_nn, CC_id, DD_id
  DOUBLE COMPLEX, DIMENSION(MPS%link_dim(pos),MPS%link_dim(pos)) :: C_id_p, C_n_p, D_id_p, D_n_p
  DOUBLE COMPLEX, DIMENSION(MPS%link_dim(pos),MPS%link_dim(pos)) :: E_id, E_n, F_id, F_n, E_temp, F_temp
  DOUBLE COMPLEX, DIMENSION(MPS%link_dim(pos),MPS%link_dim(pos)) :: LL1, LL2
  DOUBLE COMPLEX, DIMENSION(MPS%link_dim(pos),MPS%link_dim(pos)) :: Temp1, Temp2
  DOUBLE COMPLEX, DIMENSION(MPS%link_dim(pos-1),MPS%link_dim(pos+1)) :: Temp3
  DOUBLE COMPLEX, DIMENSION(MPS%link_dim(pos),MPS%link_dim(pos-1))   :: Temp4
  INTEGER, DIMENSION(MPS%link_dim(pos)) :: ipiv1, ipiv2
  
  LOGICAL Initial
  INTEGER ii, jj, kk, ww, m, d, counter, counter2
  DOUBLE PRECISION N, N_exact, lambda, norm, deriv, aux1, aux2, max_err
  
  d = MPS%local_dim
  m = MPS%link_dim(pos)
  
  norm = 0.d0
  N_exact = 0.d0
  ! Calculate norm and exact expectation value of the particle number operator with the MPS before truncating
  DO ii = 1, d
    DO jj = 1, d
      aux1 = SUM(BP(ii,jj)%M*CONJG(BP(ii,jj)%M))
      norm = norm + aux1
      N_exact = N_exact + (ii+jj-2.0d0)*aux1
    ENDDO
  ENDDO
  ! Normalize the MPS
  DO ii = 1, d
    DO jj = 1, d
      BP(ii,jj)%M = BP(ii,jj)%M/SQRT(norm)
    ENDDO
  ENDDO
  N_exact = N_exact/norm

  ! Initialize the procedure
  ! Transfer the MPS matrices to A1 and A2 for convenience
  ! A1 and A2 are the matrix called B^[n_i] and B^[n_{i+1}] in the paper
  ! 'prev' stands for previous
  DO ii = 1, d
    ALLOCATE(A1(ii)%M(MPS%link_dim(pos-1),m),A1_prev(ii)%M(MPS%link_dim(pos-1),m),M1(ii)%M(MPS%link_dim(pos-1),m))
    ALLOCATE(A2(ii)%M(m,MPS%link_dim(pos+1)),A2_prev(ii)%M(m,MPS%link_dim(pos+1)),M2(ii)%M(m,MPS%link_dim(pos+1)))
    A1_prev(ii)%M = MPS%site(pos)%A(ii)%M
    A2_prev(ii)%M = MPS%site(pos+1)%A(ii)%M 
  ENDDO

  ! Calculate the matrices C_id, D_id, C_n, D_n (see paper)
  ! C_nn and D_nn are defined as C_n and D_n respectively but with n**2 in place of n
  ! the suffix '_p' stands for previous
  C_id_p = 0.
  D_id_p = 0.
  C_n_p  = 0.
  D_n_p  = 0.
  C_nn   = 0.
  D_nn   = 0.
  DO ii = 1, d
    Temp1 = MATMUL(TRANSPOSE(CONJG(A1_prev(ii)%M)),A1_prev(ii)%M)
    C_id_p = C_id_p + Temp1
    C_n_p  = C_n_p  + (ii-1.0d0)*Temp1
    C_nn = C_nn + (ii-1.0d0)**2*Temp1
    Temp2 = MATMUL(A2_prev(ii)%M,TRANSPOSE(CONJG(A2_prev(ii)%M)))
    D_id_p = D_id_p + Temp2
    D_n_p  = D_n_p  + (ii-1.0d0)*Temp2
    D_nn = D_nn + (ii-1.0d0)**2*Temp2
  ENDDO
  norm = Trace(m,C_id_p,D_id_p)
  N = (Trace(m,C_n_p,D_id_p) + Trace(m,C_id_p,D_n_p))/norm
  
  ! If the truncated and exact particle-number operator expectation value coincide exit the subroutine
  IF (ABS(N_exact-N) <= 1.0d-12) RETURN
  ! Calculate the derivative of the particle number expectation value as a function of lambda (mu_2 in the paper)
  ! at the point lambda = 0. The formula presented in the paper simplifies for lambda = 0 and the simplified form is used below.
  aux1 = 0.d0
  aux2 = 0.d0
  DO ii = 1, m
    DO jj = 1, m
      aux1 = aux1 + C_n_p(ii,jj)*S(jj)*C_n_p(jj,ii)/S(ii)
      aux2 = aux2 + D_n_p(ii,jj)*S(jj)*D_n_p(jj,ii)/S(ii)
    ENDDO
  ENDDO
  deriv  = 2*(Trace(m,C_nn,D_id_p)+Trace(m,C_id_p,D_nn)+4*Trace(m,C_n_p,D_n_p)+aux1+aux2)/norm - 4*N**2 ! the derivative with the simplified formula is used here
  lambda = (N_exact-N)/deriv
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   print *, 'Initialization'
!   print *, 'Norm = ', norm
!   print *, 'N_exact', n_exact
!   print *, 'N_exact - N', n_exact - n
!   print *, 'deriv = ', deriv
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Start of the algorithm iterations
  counter = 1		! count the number of iterations of the algorithm
  Initial = .TRUE.	! flag that indicates that a new iteration has started
  main: DO
    ! Multiply the current matrix A1 and A2 by the time-evolved matrix before truncation called BP (M in the paper) (BP is a collection of pointers to B)
    ! note the use of the "previous" A1_prev and A2_prev. A1_prev, A2_prev are initialized above.
    DO ii = 1, d
      M1(ii)%M = 0.d0
      M2(ii)%M = 0.d0
      DO kk = 1, d
	M1(ii)%M = M1(ii)%M + MATMUL(BP(ii,kk)%M,TRANSPOSE(CONJG(A2_prev(kk)%M))) 
	M2(ii)%M = M2(ii)%M + MATMUL(TRANSPOSE(CONJG(A1_prev(kk)%M)),BP(kk,ii)%M)
      ENDDO
    ENDDO

    ! the new A1, A2 are a function of lambda and the Newton-Ralphson method is used to find the correct value of lambda
    ! which satisfies particle number conservation
    ! begin the Newton-Ralphson algorithm
    counter2 = 0
    newton: DO
      ! If this is the first iteration of the Newton-Ralphson method use the value of lambda calculated in the previous iteration
      ! of the algorithm otherwise update lambda using the derivative of the particle number expectation value
      IF (Initial) THEN 
	Initial = .FALSE.
      ELSE
	lambda = lambda + (N_exact-N)/deriv
      ENDIF
      ! Update the matrices C_id, C_n, D_id, D_n, E_id, E_n, F_id, F_n, A1, A2
      ! using the updated value of lambda
      ! Note the use of C_id_p, D_id_p, C_n_p, D_n_p
      C_id = 0.
      C_n  = 0.
      D_id = 0.
      D_n  = 0.
      E_id = 0.
      E_n  = 0.
      F_id = 0.
      F_n  = 0.
      DO ii = 1, d     
	E_temp = (ii-1.0d0)*D_id_p + D_n_p
	LL1 = D_id_p - lambda*E_temp
	CALL Factorize(m,LL1,ipiv1) ! LDL^T decomposition for calculating the inverse
	CALL Solve(m,m,LL1,ipiv1,E_temp) ! apply the inverse of LL1 to E_temp

	F_temp = (ii-1.0d0)*C_id_p + C_n_p
	LL2 = C_id_p - lambda*F_temp
	CALL Factorize(m,LL2,ipiv2) ! LDL^T decomposition for calculating the inverse
	CALL Solve(m,m,LL2,ipiv2,F_temp) ! apply the inverse of LL2 to F_temp

	Temp4 = TRANSPOSE(CONJG(M1(ii)%M))
	CALL Solve(m,MPS%link_dim(pos-1),LL1,ipiv1,Temp4) ! Apply the inverse of LL1 to M1 (M1 = M*B^{+[n_{i+1}]} in the paper)
	A1(ii)%M = TRANSPOSE(CONJG(Temp4))	! New value of A1(lambda) (\tilde B^[n_i](\mu_2) ) in the paper) calculated using the new value of lambda 

	A2(ii)%M = M2(ii)%M
	CALL Solve(m,MPS%link_dim(pos+1),LL2,ipiv2,A2(ii)%M) ! Apply the inverse of LL2 to M2 (M2 = B^{+[n_{i}]}*M in the paper)
							     ! New value of A2(lambda) (\tilde B^[n_{i+1}](\mu_2) in the paper) calculated using the new value of lambda 
	
	! Calculate C_id, D_id, C_n, D_n, E_id, E_n, F_id, F_n
	Temp1 = MATMUL(TRANSPOSE(CONJG(A1(ii)%M)),A1(ii)%M)
	C_id = C_id + Temp1
	C_n  = C_n  + (ii-1.0d0)*Temp1 
	E_temp = MATMUL(E_temp,Temp1)
	E_temp = E_temp + TRANSPOSE(CONJG(E_temp))
	E_id = E_id + E_temp
	E_n  = E_n  + (ii-1.0d0)*E_temp
	
	Temp2 = MATMUL(A2(ii)%M,TRANSPOSE(CONJG(A2(ii)%M)))
	D_id = D_id + Temp2
	D_n  = D_n  + (ii-1.0d0)*Temp2
	F_temp = MATMUL(F_temp,Temp2)
	F_temp = F_temp + TRANSPOSE(CONJG(F_temp))
	F_id = F_id + F_temp
	F_n =  F_n + (ii-1.0d0)*F_temp
      ENDDO
      norm = Trace(m,C_id,D_id)	! Norm of the MPS with the new value of lambda
      N = (Trace(m,C_n,D_id) + Trace(m,C_id,D_n))/norm	! particle number expectation value with the new value of lambda
      deriv = (Trace(m,E_n,D_id)+Trace(m,C_n,F_id)+Trace(m,E_id,D_n)+Trace(m,C_id,F_n))/norm - &
		    N*(Trace(m,E_id,D_id)+Trace(m,C_id,F_id))/norm	! derivative of the particle number expectation value calculated with the new value of lambda
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       print *, 'newton step ', counter2, '---------------------'
!       print *, 'norm =', norm 
!       print *, 'N_exact - N', N_exact - N
!       print *, 'deriv = ', deriv
!       print *, 'lambda = ', lambda
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (counter2 >= 10) THEN
	STOP 'Maximum number of iterations for Newton-Ralphson root-finding exceeded'
      ELSEIF (ABS(N_exact - N) <= 1.0d-14)  THEN ! If the expectation value of the number particle with the new value of lambda is closed to the exact one stop the Newton-Ralphson algorithm
	! normalize the various matrices 
	C_id = C_id/SQRT(norm)
	C_n  = C_n/SQRT(norm)
	D_id = D_id/SQRT(norm)
	D_n  = D_n/SQRT(norm)
	DO kk = 1, d
	  A1(kk)%M = A1(kk)%M/norm**0.25d0
	  A2(kk)%M = A2(kk)%M/norm**0.25d0
	ENDDO
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	print *, 'Norm A1-A2', Normalization(d,m,A1,A2)
!     	print *, 'C_id',  sum(ABS(C_id))
!     	print *, 'C_n', sum(ABS(C_n))
!     	print *, 'D_id',  sum(ABS(D_id))
!     	print *, 'D_n',   sum(ABS(D_n))
!     	print *, 'Trace(C_id,D_id)', REAL(Trace(m,C_id,D_id))
!     	print *, 'N1+N2', REAL(Trace(m,C_n,D_id) + Trace(m,C_id,D_n))
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	EXIT newton
      ELSE ! next iteration of Newton-Ralphson method
	counter2 = counter2 + 1
      ENDIF
    ENDDO newton
    
    ! Calculate the distance between the normalized A1,A2 matrices and the previous ones 
    Temp1 = 0.0d0
    Temp2 = 0.0d0
    DO kk = 1, d
      Temp1 = Temp1 + MATMUL( TRANSPOSE(CONJG(A1_prev(kk)%M)), A1(kk)%M )
      Temp2 = Temp2 + MATMUL( A2(kk)%M, TRANSPOSE(CONJG(A2_prev(kk)%M)) )
    ENDDO
!    print *, '------'
!    print *, ABS(Trace(m,Temp1,Temp2))
    max_err = 1.0d0-ABS(Trace(m,Temp1,Temp2))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    print *, '-----------------------------------'
!    print *, 'max_err', max_err
!    print *, '------'
!    print *, '                     '
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (max_err <= 1.0d-12 .OR. counter >= 30) THEN ! termination condition for the full correction algorithm 
      DO ii = 1, d
	MPS%site(pos)%A(ii)%M   = A1(ii)%M
	MPS%site(pos+1)%A(ii)%M = A2(ii)%M
      ENDDO
      IF (counter >= 30) print *, 'Maximum number of interations exceeded', max_err
!      CALL Test(d,m,MPS%site(pos)%A,MPS%site(pos+1)%A,norm,N)
!      print *, 'Normalization', norm 
!      print *, 'Particle number difference', ABS(N_exact-N)
!      print *, 'N_exact' ,N_exact, 'N', N
!      IF (ABS(N_exact-N) > 2.0d-14) STOP 'Error in particle number correction'  
      EXIT main
    ELSE ! Next iteration, start with the A1, A2 matrices calculated using the new value of lambda
      counter = counter + 1
      DO ii = 1, d
	A1_prev(ii)%M = A1(ii)%M
	A2_prev(ii)%M = A2(ii)%M
      ENDDO
      C_id_p = C_id
      D_id_p = D_id
      C_n_p  = C_n
      D_n_p  = D_n
      Initial = .TRUE.
    ENDIF
  ENDDO main
  
  DO ii = 1, d
    DEALLOCATE(A1(ii)%M,A1_prev(ii)%M,M1(ii)%M,A2(ii)%M,A2_prev(ii)%M,M2(ii)%M)
  ENDDO
  
  CONTAINS

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!   SUBROUTINE Test(dd,cd,AA,BB,mm,NN)
!     
!     IMPLICIT NONE
!     INTEGER, INTENT(IN) :: dd, cd
!     TYPE(Matrix), DIMENSION(dd), INTENT(IN) :: AA
!     TYPE(Matrix), DIMENSION(dd), INTENT(IN) :: BB
!     DOUBLE PRECISION, INTENT(OUT) :: mm, NN
!     DOUBLE COMPLEX, DIMENSION(cd,cd) :: Temp1, Temp2, Temp1_N, Temp2_N
!     DOUBLE COMPLEX Normalization
!     INTEGER ii
!     
!     Temp1 = 0.d0
!     Temp2 = 0.d0
!     Temp1_N = 0.d0
!     Temp2_N = 0.d0
!     DO ii = 1, dd
!       Temp1 = Temp1 + MATMUL( TRANSPOSE(CONJG(AA(ii)%M)), AA(ii)%M )
!       Temp2 = Temp2 + MATMUL( BB(ii)%M, TRANSPOSE(CONJG(BB(ii)%M)) )
!       Temp1_N = Temp1_N + (ii-1.0d0)*MATMUL( TRANSPOSE(CONJG(AA(ii)%M)), AA(ii)%M )
!       Temp2_N = Temp2_N + (ii-1.0d0)*MATMUL( BB(ii)%M, TRANSPOSE(CONJG(BB(ii)%M)) )
!     ENDDO
!     mm = Trace(cd,Temp1,Temp2)
!     NN = Trace(cd,Temp1_N,Temp2) + Trace(cd,Temp1,Temp2_N)
!     
!   END SUBROUTINE Test
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  PURE FUNCTION Trace(dd,A,B)
    ! Calculate the trace Tr[AB] of m x m matrices A and B
    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dd
    DOUBLE COMPLEX, DIMENSION(dd,dd), INTENT(IN) :: A, B
    DOUBLE COMPLEX Trace
    INTEGER ii, jj
    
    Trace = 0.d0
    DO ii = 1, dd
      DO jj = 1, dd
	Trace = Trace + A(ii,jj)*B(jj,ii)
      ENDDO
    ENDDO
  
  END FUNCTION Trace
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  SUBROUTINE Factorize(m,A,ipiv)
    ! Calculate the LDL^T decompoition of A
    ! Result is stored in A and ipiv
  
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: m
    DOUBLE COMPLEX, DIMENSION(m,m), INTENT(INOUT) :: A
    INTEGER, DIMENSION(m), INTENT(OUT) :: ipiv
    DOUBLE COMPLEX, DIMENSION(10*m) :: Work
    INTEGER info
    
    CALL ZHETRF('L',m,A,m,ipiv,Work,10*m,info)
    IF (info /= 0) STOP 'Error in Cholesky Factorization'
  
  END SUBROUTINE Factorize
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  SUBROUTINE Solve(m,nc,A,ipiv,B)
    ! Apply the inverse of A (m x m) to B (m x nc)
    ! Use A and ipiv from Factorize (LDL^T) decomposition
    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: m, nc
    DOUBLE COMPLEX, DIMENSION(m,m), INTENT(IN) :: A
    DOUBLE COMPLEX, DIMENSION(m,nc), INTENT(INOUT) :: B
    INTEGER, DIMENSION(m), INTENT(IN) :: ipiv
    INTEGER info
        
    CALL ZHETRS('L',m,nc,A,m,ipiv,B,m,info)
    IF (info /= 0) STOP 'Error in Solving linear system'
  
  END SUBROUTINE Solve
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE Particle_Number_Conservation_Correction


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Measure_Local_Observable(MPS,op,out)
  ! Measure the local operator op on the MPS and store the result in out

  IMPLICIT NONE
  TYPE(Matrix_Product_State), INTENT(INOUT) :: MPS
  TYPE(One_Site_Operator), INTENT(IN) :: op
  DOUBLE COMPLEX, DIMENSION(MPS%lattice_length), INTENT(OUT) :: out
  TYPE(Matrix), DIMENSION(MPS%lattice_length) :: Temp
  INTEGER ii, jj, r, c, l
  DOUBLE COMPLEX val
  
  CALL Normalize_MPS(MPS,'R')
  out = 0.
  l = MPS%lattice_length
  DO ii = 1, l-1
    ALLOCATE(Temp(ii)%M(MPS%link_dim(ii),MPS%link_dim(ii)))
    Temp(ii)%M = 0.
    DO jj = 1, op%num_Elements
      r   = op%element(jj)%r
      c   = op%element(jj)%c
      val = op%element(jj)%v
      Temp(ii)%M = Temp(ii)%M + val*MATMUL(TRANSPOSE(CONJG(MPS%site(ii)%A(r)%M)),MPS%site(ii)%A(c)%M)
    ENDDO
    DO jj = 1, MPS%link_dim(ii)  !Trace
      out(ii) = out(ii) + Temp(ii)%M(jj,jj)
    ENDDO
    CALL MPS_Singular_Value_Decomposition(MPS,ii,'R')
  ENDDO
  ALLOCATE(Temp(l)%M(MPS%link_dim(l),MPS%link_dim(l)))
  Temp(l)%M = 0.
  DO jj = 1, op%num_Elements
    r   = op%element(jj)%r
    c   = op%element(jj)%c
    val = op%element(jj)%v
    Temp(l)%M = Temp(l)%M + val*MATMUL(TRANSPOSE(CONJG(MPS%site(l)%A(r)%M)),MPS%site(l)%A(c)%M)
  ENDDO
  DO jj = 1, MPS%link_dim(l)  	!Trace
      out(l) = out(l) + Temp(l)%M(jj,jj)
  ENDDO
  DO ii = 1, l
    DEALLOCATE(Temp(ii)%M)
  ENDDO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Find_Blocks(d,op,binf)
  ! Routine to find the block structure of the two-site operator op
  ! The result is stored in binf

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: d 
  TYPE(Two_Site_Operator), INTENT(IN) :: op
  TYPE(Block_Info_Type), INTENT(OUT)  :: binf
  INTEGER, DIMENSION(d**2):: temp	
  INTEGER, DIMENSION(:), ALLOCATABLE :: temp_ind
  INTEGER ii, jj, counter, r1, c1, r2, c2
  LOGICAL Changed 

   binf%local_dim = d
  ALLOCATE(binf%map(d,d))
  counter = 0
  DO ii = 1, d
    DO jj = 1, d
      counter = counter + 1
      binf%map(ii,jj)%blk = counter
      binf%map(ii,jj)%ind = 0
    ENDDO
  ENDDO
  counter = 0
  Changed = .TRUE.
  DO WHILE(Changed)
    counter = counter + 1
    Changed = .FALSE.
    DO ii = 1, op%num_Elements
      r1 = op%element(ii)%r1
      c1 = op%element(ii)%c1
      r2 = op%element(ii)%r2
      c2 = op%element(ii)%c2
      IF (binf%map(r1,r2)%blk /= binf%map(c1,c2)%blk) THEN
	binf%map(r1,r2)%blk = MIN(binf%map(r1,r2)%blk,binf%map(c1,c2)%blk)
	binf%map(c1,c2)%blk = MIN(binf%map(r1,r2)%blk,binf%map(c1,c2)%blk)
	Changed = .TRUE.
      ENDIF
    ENDDO
  ENDDO
  temp = 0
  DO ii = 1, d
    DO jj = 1, d
      temp(binf%map(ii,jj)%blk) = temp(binf%map(ii,jj)%blk) + 1
    ENDDO
  ENDDO
  counter = 0
  DO ii = 1, d**2
    IF (temp(ii) /= 0) counter = counter + 1
  ENDDO
  binf%num_blk = counter
  ALLOCATE(binf%blk_dim(binf%num_blk))
  counter = 0
  DO ii = 1, d**2
    IF (temp(ii) /= 0) THEN 
      counter = counter + 1
      binf%blk_dim(counter) = temp(ii)
      temp(ii) = counter
    ENDIF
  ENDDO
  DO ii = 1, d
    DO jj = 1, d
      binf%map(ii,jj)%blk = temp(binf%map(ii,jj)%blk)
    ENDDO  
  ENDDO
  ALLOCATE(temp_ind(binf%num_blk))
  temp_ind = 0
  DO ii = 1, d
    DO jj = 1, d
      temp_ind(binf%map(ii,jj)%blk) = temp_ind(binf%map(ii,jj)%blk) + 1
      binf%map(ii,jj)%ind = temp_ind(binf%map(ii,jj)%blk)
    ENDDO
  ENDDO
  ALLOCATE(binf%block(binf%num_blk))
  DO ii = 1, binf%num_blk
    ALLOCATE(binf%block(ii)%v(binf%blk_dim(ii)))
  ENDDO
  DO ii = 1, d
    DO jj = 1, d
      binf%block(binf%map(ii,jj)%blk)%v(binf%map(ii,jj)%ind)%ind1 = ii
      binf%block(binf%map(ii,jj)%blk)%v(binf%map(ii,jj)%ind)%ind2 = jj
     ENDDO
  ENDDO
  DEALLOCATE(temp_ind)
!   11 FORMAT('('I4,I4,')   (',I4,I4,')')
!   DO ii = 1, d
!     DO jj = 1, d
!       print 11, ii, jj, binf%map(ii,jj)%blk, binf%map(ii,jj)%ind
!     ENDDO
!   ENDDO 
!   print *, '------------------------------------------------------------------'
!   DO ii = 1, binf%num_blk
!     DO jj = 1, binf%blk_dim(ii)
!       print 11, binf%block(ii)%v(jj)%ind1, binf%block(ii)%v(jj)%ind2, ii, jj
!     ENDDO
!   ENDDO 
  
END SUBROUTINE Find_Blocks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Exponentiate(op_in,binf,coeff,op_out)
  ! Exponentiate the two-site operator op_in into the two_site operator op_out
  ! the coefficient of the exponent is coeff, and binf provides informations on the block structure in order to exponentiate 
  ! block by block
  
  IMPLICIT NONE
  TYPE(Two_Site_Operator), INTENT(IN) :: op_in
  TYPE(Two_Site_Operator), INTENT(OUT) :: op_out
  TYPE(Block_Info_Type), INTENT(IN) :: binf
  DOUBLE COMPLEX, INTENT(IN) :: coeff
  TYPE(Matrix), DIMENSION(:), ALLOCATABLE :: blk
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: W
  DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: Temp
  INTEGER ii, jj, kk, r1, c1, r2, c2, counter
  
  ALLOCATE(blk(binf%num_blk))
  DO ii = 1, binf%num_blk
    ALLOCATE(blk(ii)%M(binf%blk_dim(ii),binf%blk_dim(ii)))
    blk(ii)%M = 0.
  ENDDO
  DO ii = 1, op_in%num_Elements
    r1  = op_in%element(ii)%r1
    c1  = op_in%element(ii)%c1
    r2  = op_in%element(ii)%r2
    c2  = op_in%element(ii)%c2
    IF (binf%map(r1,r2)%blk /= binf%map(c1,c2)%blk) STOP 'Error in Exponentiate'
    blk(binf%map(c1,c2)%blk)%M(binf%map(r1,r2)%ind,binf%map(c1,c2)%ind) = op_in%element(ii)%v
  ENDDO
  DO ii = 1, binf%num_blk
    ALLOCATE(W(binf%blk_dim(ii)),Temp(binf%blk_dim(ii),binf%blk_dim(ii)))
    CALL DiagMatrix(binf%blk_dim(ii),blk(ii)%M,W)
    DO jj = 1, binf%blk_dim(ii)
      Temp(:,jj) = blk(ii)%M(:,jj)*EXP(coeff*W(jj))
    ENDDO
    blk(ii)%M = MATMUL(Temp,TRANSPOSE(CONJG(blk(ii)%M)))
    DEALLOCATE(Temp,W)
  ENDDO
  counter = 0
  DO ii = 1, binf%num_blk
    DO jj = 1, binf%blk_dim(ii)
      DO kk = 1, binf%blk_dim(ii)
	IF (ABS(blk(ii)%M(jj,kk)) >= 1.0d-12) counter = counter + 1
      ENDDO
    ENDDO
  ENDDO
  op_out%num_Elements = counter
  ALLOCATE(op_out%element(op_out%num_Elements))
  counter = 0.
  DO ii = 1, binf%num_blk
    DO jj = 1, binf%blk_dim(ii)
      DO kk = 1, binf%blk_dim(ii)
	IF (ABS(blk(ii)%M(jj,kk)) >= 1.0d-12) THEN
	  counter = counter + 1
	  op_out%element(counter)%r1 = binf%block(ii)%v(jj)%ind1
	  op_out%element(counter)%c1 = binf%block(ii)%v(kk)%ind1
	  op_out%element(counter)%r2 = binf%block(ii)%v(jj)%ind2
	  op_out%element(counter)%c2 = binf%block(ii)%v(kk)%ind2
	  op_out%element(counter)%v = blk(ii)%M(jj,kk)
	ENDIF
      ENDDO
    ENDDO
    DEALLOCATE(blk(ii)%M)
  ENDDO
  IF (counter /= op_out%num_Elements) STOP 'Error in Exponentiate'
  DEALLOCATE(blk)
  
  CONTAINS
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  SUBROUTINE DiagMatrix(n,M,W) ! Uses lapack routines to diagonalize Hermitian matrices

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE COMPLEX, DIMENSION(n,n), INTENT(IN) :: M
  DOUBLE PRECISION, INTENT(OUT), DIMENSION(n) :: W
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

END SUBROUTINE Exponentiate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
END MODULE Open_DMRG