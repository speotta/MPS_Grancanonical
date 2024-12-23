!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!   Modified TDMRG algorithm for grand-canonical MPS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!	This is a Fortran 95 implementation of a Density Matrix Renormalization Group Algorithm using a Matrix Product State approximation of the wavefunction.
!!!
!!!	The MPS is optimized and evolved in time in the grandcanonical ensemble, which means that the MPS is not invariant under the Abelian symmetry 
!!!	generated by the total number of particles.
!!!
!!!	However a modified version of the usual TDMRG truncation procedure helps to ensure the conservation of the expectation value 
!!!	of the number of particles as described in the paper "Improving the Gutzwiller Ansatz with Matrix Product States" 
!!!	arXiv:1307.8416 by S. Peotta and M. Di Ventra. If you use and modify the code for you work, please provide a reference to the above paper.
!!!
!!!	The code is not documented so far (April 2014) but it commented thorughout and should be easy to understand if you have knowledge of MPS-based algorithm. 
!!!	The code can simulate only the Bose-Hubbard Hamiltonian with arbitrary potentials, but it is easy to implement arbitrary models 
!!!	with nearest-neighbour interactions. 
!!!
!!!	The code has been developed an tested on the Ubuntu 14.04 (Trusty Tahr) operating system and requires the GNU GFortran compiler and 
!!!	the standard LAPACK and BLAS library or, alternatively, the ifort compiler with Math Kernel Library in place of Lapack and Blas
!!!
!!!	To compile the code simply run the makefile using the command
!!!	$ make
!!!
!!!	to clean up use the command
!!!	$ make clean 
!!!
!!!	If you have questions or you would like to provide your own improvements to the code 
!!!	feel free to contact me (see my homepage http://physics.ucsd.edu/~speotta/).
!!!
!!!	22th April 2014
!!!
!!!	Sebastiano Peotta, Post-Doc in Condensed Matter Physics at University of California-San Diego
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM mps_grand

  USE Bose_Hubbard

  IMPLICIT NONE

  ! Problem parameters
  INTEGER l			! lattice length
  INTEGER d			! dimension of the site Hilbert space
  INTEGER steps			! number of steps of real time 
  INTEGER osteps		! number of steps between two evalutations of the observables
  INTEGER sweeps_White_Correction	! number of sweeps of single site optimization with White correction
  INTEGER sweeps_Double_Site		! number of sweeps of double site optimization
  INTEGER sweeps_Single_Site		! number of sweeps of single site optimization
  INTEGER m_Min				! minimum value of the link dimension
  INTEGER m_Max				! maximum value fo the link dimension
  DOUBLE PRECISION u			! Hubbard interaction
  DOUBLE PRECISION mu 			! chemical potential
  DOUBLE PRECISION phi			! phase difference in the phase quench
  DOUBLE PRECISION dt			! time step
  DOUBLE PRECISION eps			! truncation error
  DOUBLE PRECISION ground_state_erg	! ground state energy
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: V0, Vt	! initial potential (V0) and potential for time evolution (Vt)
  
  TYPE(Matrix_Product_State) :: MPS			! Matrix Product State
  DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE :: out	! output vector with the local density
  TYPE(One_Site_Operator) :: n				! number operator
  TYPE(Hamiltonian) :: time_Ham				! Hamiltonian

  ! working local variables
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: psi
  INTEGER ii
  REAL time1, time2, time3, time4
  
  ! Read the parameters from an external file
  OPEN(unit=100, file='./mps_grand.in', status = 'OLD')
  READ(100,*) l
  READ(100,*) d
  READ(100,*) u
  READ(100,*) mu
  READ(100,*) phi
  READ(100,*) dt
  READ(100,*) steps
  READ(100,*) osteps  
  READ(100,*) eps
  READ(100,*) sweeps_White_Correction
  READ(100,*) sweeps_Double_Site
  READ(100,*) sweeps_Single_Site
  READ(100,*) m_Min
  READ(100,*) m_Max
  CLOSE(unit=100)
  
  IF (MOD(l,2) == 1) STOP 'Lattice length must be even'
  
  ! Read the external potential for the static optimization
  ALLOCATE(V0(l))
  OPEN(unit=101, file='./potential_init.in', status = 'OLD')
  READ(101,*) V0
  CLOSE(unit = 101)
 
  ! Read the external potential for the time evolution
  ALLOCATE(Vt(l))
  Vt = 0.d0
  OPEN(unit=102, file='./potential_time.in', status = 'OLD')
  READ(102,*) Vt
  CLOSE(unit = 102)
 
  ! Output the parameters on the standard output
  WRITE (*,10) l, d, u, mu, phi, dt, steps, osteps, eps, sweeps_White_Correction, &
		sweeps_Double_Site, sweeps_Single_Site, m_Min, m_Max
  10 FORMAT (/,T20,'Number conserving grancanonical MPS',//,//,'Parameters',//,  &
  1X, 'L    =', I6,  /  & 
  1X, 'D    =', I6,  /  &
  1X, 'U    =', ES10.3,/  &
  1X, 'MU   =', ES10.3,/  &
  1X, 'Phi  =', ES10.3,/  &
  1X, 'dt   =', ES10.3,/  &
  1X, 'stp  =', I6,  /  &
  1X, 'ostp =', I6,  /  &
  1X, 'eps  =', ES10.3,/ &
  1X, 'sWC  =', I6,  /, &
  1X, 'sDS  =', I6,  /, &
  1X, 'sSS  =', I6,  /, &
  1X, 'mMin =', I6,  /, &
  1X, 'mMax =', I6,  /)

  ! Provide a translational invariant Gutzwiller wavefunction as an initial guess
  ! using the Bose-Hubbard Hamiltonian with the given values of the interaction (u)
  ! and chemical potential (mu)
  ! The result is stored in psi(:)
  ALLOCATE(psi(d))
  CALL Gutzwiller_Translational_Invariant_Guess(u,mu,d,psi)

  ! Inizialize a MPS wavefunction with the traslational invariant Guztwiller wavefunction
  CALL Initialize_MPS_Bose_Hubbard(d,psi,l,MPS)
  
  ! Alternatively use a random guess
  ! CALL Initialize_MPS_Random(d,l,1,MPS)

  ! Build the Bose-Hubbard Hamiltonian (Ham) needed to find the ground state
  CALL Build_Bose_Hubbard_Hamiltonian(l,d,u,mu,0.d0,V0,Ham)
  
  ! Initialize the effective Hamiltonian (CALL CPU_time(time3)Eff_Ham) using the current MPS and the Hamiltonian (Ham)
  CALL Initialize_Effective_Hamiltonian(MPS,Ham,Eff_Ham)
  
  ! Begin the optimization in order to find the ground state
  ! Open the file where the expectation value of the energy is stored during the sweeping process
  ! useful to check that the energy is decreasing as the variational wavefunction is optimized
  OPEN(unit=104, file='./energy.dat', status = 'REPLACE')
  ! Single site optimization with White correction, number of sweeps = sweeps_White_Correction
  CALL CPU_time(time1)
  CALL MPS_Jacobi_Davidson_Optimization_White_Correction(MPS,sweeps_White_Correction,m_Min,m_Max,eps)
  ! Single site optimization with White correction, number of sweeps = sweeps_Double_Site 
  CALL MPS_Jacobi_Davidson_Optimization_Double_Site(MPS,sweeps_Double_Site,m_Min,m_Max,eps)
  ! Single site optimization, number of sweeps = sweeps_Single_Site
  CALL MPS_Jacobi_Davidson_Optimization_Single_Site(MPS,sweeps_Single_Site)
  CALL CPU_time(time2)
  print *, 'Time for ground state optimization, t = ', time2-time1 
  CLOSE(unit=104)
  
  ! Build the number operator for a single site
  CALL Build_Number_Operator(d,n)
  ! Measure the density, the result is stored in out
  ALLOCATE(out(l))
  CALL Measure_Local_Observable(MPS,n,out)

  ! Write the measured value of the density in the output file
  20 FORMAT(1000ES10.3)
  OPEN(unit=103, file='./density.dat', status = 'REPLACE')
  WRITE(103,20) DREAL(out)
  CLOSE(unit = 103)
  
  ! Nprmalize the MPS
  CALL Normalize_MPS(MPS,'R')
  
  IF (steps > 0) THEN
    ! Build the Hamiltonian for the real time evolution
    CALL Build_Bose_Hubbard_Hamiltonian_Time_Evolution(l,d,u,mu,phi,Vt,time_Ham)
    ! Perform the real time evolution
    CALL CPU_time(time3)
    CALL TimeEvolution(MPS,time_Ham,.FALSE.,6,dt,steps,osteps,m_Min,m_Max,eps)
    CALL CPU_time(time4)
    print *, 'Time for real-time evolution, t = ', time4-time3
  ENDIF
  CALL CPU_time(time2)
  print *, 'Total time elapsed, t = ', time2-time1
  
END PROGRAM mps_grand

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
