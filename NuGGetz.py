#!/usr/bin/env python

"""NuGGetz is the implementation of a simple idea to correct the time evolution of a 
grancanonical Gutzwiller ansatz in order to enforce the conservation of the expectation 
value of the particle number operator. The algorithm is a modification of the well-known 
Time-Dependent Density Matrix Renormalization Group (TDMRG) algorithm and it can be 
readily generalized to a more general Matrix Product State (MPS) that does not explicitly 
conserve the particle number or other operators that generate abelian symmetries of the 
system Hamiltonian. For more details on the algorithm see S. Peotta and M. Di Ventra, 
to be published. The code, provided here as a Python module, is a chep alternative to the 
full TDMRG algorithm based on explicitly number-conserving MPSs which requires a 
substantial programming effort. We checked NuGGetz simulations against quasi-exact 
TDMRG simulations based on number conserving MPS and in all cases we found that the 
gross features of the dynamics are well captured by NuGGetz.


Conditions of use:

The NuGGetz code is free for academic use only. Any published work which utilizes 
this code should cite the following reference:

S. Peotta and M. Di Ventra
Improving the Gutzwiller ansatz with Matrix Product States
arXiv:....."""

import sys
import os
from numpy import *
from numpy.linalg import eigh, svd, solve, norm
from numpy.random import rand
from scipy.optimize import brentq
import cPickle
import time


def f (ll, v1, v2, N):
	"""
		Function whose zero is the value of the Langrange multiplier \lambda^*
		that ensures conservation of the number of particles.
	"""
	dim = v1.shape[0]
	w1 = v1/(1.-ll*arange(dim))
	w1 = w1/sqrt(vdot(w1,w1))
	w2 = v2/(1.-ll*arange(dim))
	w2 = w2/sqrt(vdot(w2,w2))
	return abs(sum(arange(dim)*w1*w1.conj())+sum(arange(dim)*w2*w2.conj())) - N


def density(c):
	"""
		Returns the density calculated from the Gutzwiller wavefunction.
	"""
	n = ndarray(c.shape[0])
	for ii in range(c.shape[0]):
		n[ii] = sum(arange(c.shape[1])*abs(c[ii,:])**2)
	return n


def apply_two_sites(U,c1,c2,correction=False):
	"""
		Apply a two-site unitary operator to the wavefunction and replaces the old
		wavefunction with the states corresponding to the largest singular value. If 
		correction=True apply the correction needed to conserve the number of particles.
		For details on the algorithm see S.Peotta and M.Di Ventra...
	"""

	dim = c1.shape[0]
	n = arange(dim)
	Psi = tensordot(U,outer(c1,c2),2)
	A,S,B = svd(Psi)
	psi1 = A[:,0]
	psi2 = B[0,:]
	
	if correction:
		atol = 1.0e-12
		rtol = 1.0e-9
		N = abs(sum(n*diagonal(dot(Psi.conj(),Psi.T)))+sum(n*diagonal(dot(Psi.conj().T,Psi))))
		meanN_1 = abs(sum(n*psi1*psi1.conj()))
		meanN_2 = abs(sum(n*psi2*psi2.conj()))
		meanN2  = abs(sum((n**2)*psi1*psi1.conj()+(n**2)*psi2*psi2.conj()))
		mm = 2*(meanN2-meanN_1**2-meanN_2**2)
		guess = -(meanN_1+meanN_2-N)/mm
		res = brentq(f, 0., 2*guess, args=(psi1,psi2,N))
		psi1_prec = psi1
		psi2_prec = psi2
		psi1 = psi1/(1-res*n)
		psi1 = psi1/sqrt(vdot(psi1,psi1))
		psi2 = psi2/(1-res*n)
		psi2 = psi2/sqrt(vdot(psi2,psi2))
		count = 0
		while not (allclose(psi1,psi1_prec,atol=atol,rtol=rtol) and allclose(psi2,psi2_prec,atol=atol,rtol=rtol)):
			aux1 = dot(Psi,psi2.conj())
			aux1 = aux1/sqrt(vdot(aux1,aux1))
			aux2 = dot(psi1.conj(),Psi)				
			aux2 = aux2/sqrt(vdot(aux2,aux2))
			meanN_1 = abs(sum(n*aux1*aux1.conj()))
			meanN_2 = abs(sum(n*aux2*aux2.conj()))
			meanN2 = abs(sum((n**2)*aux1*aux1.conj()+(n**2)*aux2*aux2.conj()))
			mm = 2*(meanN2-meanN_1**2-meanN_2**2)
 			guess = -(meanN_1+meanN_2-N)/mm
 			res = brentq(f, 0., 2.*guess, args=(aux1,aux2,N))
			psi1_prec = psi1
			psi2_prec = psi2
			psi1 = aux1/(1-res*n)
			psi1 = psi1/sqrt(vdot(psi1,psi1))		
			psi2 = aux2/(1-res*n)
			psi2 = psi2/sqrt(vdot(psi2,psi2))
			count += 1
			if count > 20:
				print 'Correction not converged'
				break
	
	return psi1, psi2


def expHam(H,tau):
	"""
		Calculate the exponential of a two-site operator that is part of the Hamiltonian
	"""
	d = H.shape[0]
	aux = ndarray((d**2,d**2),dtype = complex)
	for irow in xrange(d):
		for jrow in xrange(d):
			for icol in xrange(d):
				for jcol in xrange(d):
					aux[irow*d+jrow,icol*d+jcol] = H[irow,jrow,icol,jcol]
	l, U = eigh(aux)
	l = exp(-tau*l)
	aux = dot(U,dot(diag(l),U.conj().T))
	G = ndarray(H.shape,dtype = complex)
	for irow in xrange(d):
		for jrow in xrange(d):
			for icol in xrange(d):
				for jcol in xrange(d):
					G[irow,jrow,icol,jcol] = aux[irow*d+jrow,icol*d+jcol] 
	return G


def computeEvol(H,dt,TrotterOrder):
	"""
		Compute the two site operators used for the Trotter expansion
	"""
	if TrotterOrder == 2:
		order = [0.5,1.0,0.5]
	elif TrotterOrder == 4:
		order = [0.67560359597982889,
				1.3512071919596578,
				-0.17560359597982883,
				-1.7024143839193153,
				-0.17560359597982883,
				1.3512071919596578,
				0.67560359597982889]
	elif TrotterOrder == 6:
		order = [0.39225680523877998, 
				0.78451361047755996, 
				0.51004341191845848, 
				0.23557321335935699,
				-0.47105338540975655, 
				-1.1776799841788701, 
				0.068753168252518093, 
				1.3151863206839063,
				0.068753168252518093,
				-1.1776799841788701,
				-0.47105338540975655, 
				0.23557321335935699,
				0.51004341191845848, 
				0.78451361047755996, 
				0.39225680523877998]			
	evolU = []
	for jj in range(len(order)):
		evolU.append(expHam(H,dt*order[jj]))
	return evolU	


def Optimize(psi0,MU,U):
	"""
		Calculate the initial state by optimizing the grancanonical Gutzwiller ansatz with
		respect to a given Hamiltonian
	"""

	L, dim = psi0.shape
	
	b = zeros((dim,dim))
	for ii in range(1,dim):
		b[ii-1,ii] = sqrt(ii)
	bd = zeros((dim,dim))
	for ii in range(1,dim):
		bd[ii,ii-1] = sqrt(ii)
	n = dot(bd,b)
	n2 = dot(dot(bd,n),b)

	psi0_prev = zeros((L,dim))

	ii = 0
	while not allclose(psi0,psi0_prev,atol=1.0e-14,rtol=1.0e-10):
		ii += 1
		psi0_prev = psi0.copy()
		if ii % 2 == 1:
			beta = vdot(psi0[1,:],dot(b,psi0[1,:]))
			M = 0.5*U*n2-MU*n-beta*bd-conj(beta)*b
			erg, M = eigh(M)
			psi0[0,:] = M[:,0]	
			for site in range(1,L-1):
				beta = vdot(psi0[site-1,:],dot(b,psi0[site-1,:]))+vdot(psi0[site+1,:],dot(b,psi0[site+1,:]))
				M = 0.5*U*n2-MU*n-beta*bd-conj(beta)*b
				erg, M = eigh(M)
				psi0[site,:] = M[:,0]
			beta = vdot(psi0[L-2,:],dot(b,psi0[L-2,:]))
			M = 0.5*U*n2-MU*n-beta*bd-conj(beta)*b
			erg, M = eigh(M)
			psi0[L-1,:] = M[:,0]				
		elif ii % 2 == 0:
			beta = vdot(psi0[L-2,:],dot(b,psi0[L-2,:]))
			M = 0.5*U*n2-MU*n-beta*bd-conj(beta)*b
			erg, M = eigh(M)
			psi0[L-1,:] = M[:,0]
			for site in range(L-2,0,-1):
				beta = vdot(psi0[site-1,:],dot(b,psi0[site-1,:]))+vdot(psi0[site+1,:],dot(b,psi0[site+1,:]))
				M = 0.5*U*n2-MU*n-beta*bd-conj(beta)*b
				erg, M = eigh(M)
				psi0[site,:] = M[:,0]
			beta = vdot(psi0[1,:],dot(b,psi0[1,:]))
			M = 0.5*U*n2-MU*n-beta*bd-conj(beta)*b
			erg, M = eigh(M)
			psi0[0,:] = M[:,0]

	return psi0


def twoten(op1,op2):
	""" 
		Calculate the tensor product of two operators
	"""
	d = op1.shape[0]
	aux = ndarray((d,d,d,d))
	for irow in range(d):
		for jrow in range(d):
			for icol in range(d):
				for jcol in range(d):
					aux[irow,jrow,icol,jcol] = op1[irow,icol]*op2[jrow,jcol]
	return aux


def Time_Evolution(file_name,psi,MU,U,phi,dt,steps,osteps,TrotterOrder=6,corr=False):
	""" 
		Routine for the time evolution, for the input parameters see the main routine
		nuggetz
		TrotterOrder: order of the Trotter expansion used
	"""
 
	L, dim = psi.shape
	I = eye(dim)
	b = zeros((dim,dim))
	for ii in range(1,dim):
		b[ii-1,ii] = sqrt(ii)
	bd = zeros((dim,dim))
	for ii in range(1,dim):
		bd[ii,ii-1] = sqrt(ii)
	n = dot(bd,b)
	n2 = dot(dot(bd,n),b)

	H = 0.25*U*(twoten(n2,I)+twoten(I,n2)) - 0.5*MU*(twoten(n,I)+twoten(I,n))\
		- exp(-1.0j*phi)*twoten(bd,b) - exp(1.0j*phi)*twoten(b,bd)

	H_left = 0.5*U*twoten(n2,I)+0.25*U*twoten(I,n2)\
		-MU*twoten(n,I)-0.5*MU*twoten(I,n)\
		-exp(-1.0j*phi)*twoten(bd,b)-exp(1.0j*phi)*twoten(b,bd)

	H_right = 0.25*U*twoten(n2,I)+0.5*U*twoten(I,n2)\
		-0.5*MU*twoten(n,I)-MU*twoten(I,n)\
		-exp(-1.0j*phi)*twoten(bd,b)-exp(1.0j*phi)*twoten(b,bd)

	evolU = computeEvol(H,1.0j*dt,TrotterOrder)
	evolU_left = computeEvol(H_left,1.0j*dt,TrotterOrder)
	evolU_right = computeEvol(H_right,1.0j*dt,TrotterOrder)
	
	PSI = [density(psi.copy()),]
	f = open(file_name,'w')
	cPickle.dump(PSI,f)
	f.close()
	for ii in range(1,steps+1):
 		for jj in range(len(evolU)):
 			if jj % 2 == 0:
				psi[0,:], psi[1,:] = apply_two_sites(evolU_left[jj],psi[0,:],psi[1,:],correction=corr)
 				for site in range(2,L-2,2):
 					psi[site,:], psi[site+1,:] = apply_two_sites(evolU[jj],psi[site,:],psi[site+1,:],correction=corr)
				psi[L-2,:], psi[L-1,:] = apply_two_sites(evolU_right[jj],psi[L-2,:],psi[L-1,:],correction=corr)
 			elif jj % 2 == 1:
 				for site in range(1,L-1,2):
 					psi[site,:], psi[site+1,:] = apply_two_sites(evolU[jj],psi[site,:],psi[site+1,:],correction=corr)		
		if ii % osteps == 0:
			PSI.append(density(psi.copy()))
			f = open(file_name,'w')
			cPickle.dump(PSI,f)
			f.close()


 		
	return PSI


def nuggetz(L,MU,U,phi,dt,steps,osteps,directory,corr=False):
	""" 
		Input parameters:
		L	: lattice length
		MU  : chemical potential (to fix the number of particles
		U	: Hubbard interaction
		phi : complex hopping phase
		dt 	: time step (typically U*t <= 0.1)
		steps 	: numebr of time steps
		osteps 	: number of time steps between two 'write observable to file' operations
		directory 	: directory where to save the output (here only the on-site density is implemented
		coor : boolean, if True use the correction to keep the expectation value of number operator constant,
				which is the innovative feature of NuGGetz.
	"""

 	start = time.time()
 	
	name = 'gutzwiller_L%04d_U%0.4e'
	U = 2.*U

	dim = 7
	I = eye(dim)
	b = zeros((dim,dim))
	for ii in range(1,dim):
		b[ii-1,ii] = sqrt(ii)
	bd = zeros((dim,dim))
	for ii in range(1,dim):
		bd[ii,ii-1] = sqrt(ii)
	n = dot(bd,b)
	n2 = dot(dot(bd,n),b)
	
	c = ones(dim)/sqrt(dim)	
 	for ii in range(10000): 
 		cp = c.copy()
 		if ii == 0:
 			beta = 1.
 		else:
	 		beta = vdot(c,dot(b,c))
 		M = 0.5*U*n2-MU*n-2*beta*(bd+b)
 		erg, M = eigh(M)
 		c = M[:,0]
 		if allclose(c,cp,atol=1.0e-12,rtol=1.0e-09):
 			break

	psi0 = ndarray((L,dim),dtype=complex)
	for ii in range(psi0.shape[0]):
		psi0[ii,:] = c
	psi = Optimize(psi0,MU,U)
	file_name = './'+directory+'/'+(name % (L,U/2.))+'.dat'
	psi  = Time_Evolution(file_name,psi,MU,U,phi=phi,dt=dt,steps=steps,osteps=osteps,corr=corr)

	print "Elapsed Time: %s" % (time.time() - start)


if __name__ == "__main__":
	L = int(sys.argv[1])
	MU = float(sys.argv[2])
	U = float(sys.argv[3])
	steps = int(sys.argv[4])
	dt = float(sys.argv[5])
	osteps = int(sys.argv[6])
	directory = sys.argv[7]
	corr = int(sys.argv[8])
	phi = 0.05
	print 'Starting run with U = ', U
	nuggetz(L,MU,U,phi,dt,steps,osteps,directory,corr=corr)	
	print 'Job done with parameters (U = %+0.5e)' % U
