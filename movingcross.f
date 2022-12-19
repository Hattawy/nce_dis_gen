 	subroutine MovingCross_new (P, Px, Py, Pz, MStar, E,
     *		 Q2, Eei, Hel, Polx, Poly, Polz,  NuNRS, Phi, DSig, Nu)
C
C	Subroutine to calculate the inelastic cross section of inelastic
C	electron scattering on a moving nucleon with off-shell mass
C	MStar and momentum P. Beam energy Eei, Q^2 and Hel = product of
C	nucleon polarization and beam helicity are given in the lab. 
C	Energy transfer NuNRS is given in the Nucleon Rest System. The
C	routine returns the number of counts per total Luminosity and
C	per Delta-Q^2. To use the result, one has to bin the answers
C	in Nu and divide by the number of Monte Carlo trials.
C	Requires subroutine CrossSection.
C
C	Author: S.Kuhn  25-Feb-1993
C	Updated ("_new") S.E. Kuhn 28-Jul-2009 for proper treatment of polarization
C
	implicit none
        real*4 P        ! Total momentum of generated nucleon
        real*4 Px, Py, Pz       ! Cartesian Coordinates of P
	real*4 Polx, Poly, Polz, Pol1(3), Pol2(3), Polx2, Polz2 ! Polarization vector
	real*4 MStar
	real*4 Q2, Eei, Hel, NuNRS
	real*4 DSig, Nu ! Returned by subroutine

	real*4 Wval, TanSqTh, CotRot, E
	real*4 PeiNRS(4), EeiNRS
	real*4 LorMat(4,4)
	real*4 RotMat(4,4)
        real*4 gamma, beta, betahat(3), gammone
	real*4 qNRS(4), qrotNRS(4), qvec, costhq, sinthq
	real*4 Phi,Pi,MassN,AlPiH, Pperp
	real*4 ran1
	integer Seed
	integer I,J,K
        real*4 ppart(4,5)     ! Used to calculate final electron

	common /Constants/ Pi,MassN,AlPiH
	common /Particles/ ppart
	common /randomizer/ Seed
	
	Pol1(1) = Polx
	Pol1(2) = Poly
	Pol1(3) = Polz
	dsig = 0.0
	Nu = 0.0
	do I = 1,4
	  ppart(1,I) = 0.0
	enddo
	gamma = E/MStar
	gammone = gamma - 1.0
	beta = P/E
	if (P .gt. 0.001) then
	  betahat(1) = Px/P
	  betahat(2) = Py/P
	  betahat(3) = Pz/P
	else
	  do I = 1 , 3
	    betahat(I) = 0.0
	  enddo
	endif
	LorMat(4,4) = gamma
	do I = 1 , 3
	  LorMat(4,I) = -gamma*beta*betahat(I)
	  LorMat(I,4) = LorMat(4,I)
	  do J = 1 , 3
	    LorMat(I,J) = gammone*betahat(I)*betahat(J)
	    if (I .eq. J) LorMat(I,J) = LorMat(I,J) + 1.0
	  enddo
	enddo

	do I = 1 , 4
	  PeiNRS(I) = (LorMat(I,3) + LorMat(I,4))*Eei
	enddo
	EeiNRS = PeiNRS(4)
	if((EeiNRS .le. NuNRS).or. (Q2.le.0.0))return
	if(Q2.ge.(4*EeiNRS*(EeiNRS-NuNRS)))return
	
	Wval = sqrt(MStar*(MStar+2.0*NuNRS) - Q2)
	if(Wval .lt. 0.934)return
	
	Pperp = sqrt(PeiNRS(1)*PeiNRS(1) + PeiNRS(2)*PeiNRS(2))
	if (Pperp .gt. 0.0) then
	  CotRot = PeiNRS(3)/Pperp
	  do I = 1 , 3
	    RotMat(I,1) = PeiNRS(I)/EeiNRS*CotRot
	    RotMat(I,3) = PeiNRS(I)/EeiNRS
	  enddo
          RotMat(3,1) = - Pperp/EeiNRS
	  RotMat(1,2) = - PeiNRS(2)/Pperp
	  RotMat(2,2) = PeiNRS(1)/Pperp
	  RotMat(3,2) = 0.0
	  do I = 1 , 3
	  Pol2(I) = 0.0
	    do J = 1 , 3
	      Pol2(I) = Pol2(I) + RotMat(J,I)*Pol1(J) ! Inverse of rot. matrix
	    enddo
	  enddo
	else
	  do I = 1 , 3
	    Pol2(I) = Pol1(I)
	  enddo
	endif
C	So far, polarization is relative to beam. New q-vector will point in the direction th_q, phi
C	Now we rotate it into the system where q is along the z - direction (first q_perp along -x)

	qvec = sqrt(Q2 + NuNRS*NuNRS)
	costhq = (Q2 + 2.0*EeiNRS*NuNRS)/(2.0*EeiNRS*qvec)
	if (abs(costhq) .gt. 1.0) then
	  if (abs(costhq) .le. 1.0001) then
	    costhq = costhq/(abs(costhq)+1.0e-6)
	  else
	    write(6,*) ' Error costhq: ',EeiNRS, NuNRS, Q2, qvec, costhq
	    return
	  endif
	endif
	sinthq = sqrt(1.0 - costhq*costhq)
	Polx2 = -costhq*(cos(Phi)*Pol2(1) + sin(Phi)*Pol2(2)) + sinthq*Pol2(3)
	Polz2 = sinthq*(cos(Phi)*Pol2(1) + sin(Phi)*Pol2(2)) + costhq*Pol2(3)
C!	
c	write(3,333) P, Px, Py, Pz, MStar, E, Pol1
c333	format(' called mc with: ',9f12.5)
c	write(3,333) Q2, Eei, Hel, NuNRS, Phi
C!	
C!
c	write(3,334) Wval, PeiNRS
c334	format(' W, el NRS 4-vector: ',5f12.5)
C!	
C!
c	write(3,335) Pol2
c335	format(' 1st rotated Pol vector: ',3f12.5)
C!
C!
c	write(3,336) Polx2, Polz2
c336	format(' 2nd rotated Pol vector: ',2f12.5)
C!

	call CrossSection_new(MStar, Q2, EeiNRS, Hel, DSig, NuNRS, TanSqTh,Polx2,Polz2)
	if ((Dsig .le. 0.0).or.(TanSqTh.le.0.0)) then
	  Dsig = 0.0
	  Nu = NuNRS
	  return
	endif
csk01	if (Wval .lt. 1.08) then
csk01	  DSig = DSig * MStar/0.00939 ! Correct elastic cross sections
csk01	endif
C       DSig = DSig*(1-beta*betahat(3))
C	I'm not REALLY sure what the correct form is here...
	qNRS(4) = NuNRS
	qNRS(1) = qvec*sinthq
	qNRS(2) = qNRS(1)*sin(Phi)
	qNRS(1) = qNRS(1)*cos(Phi)
	qNRS(3) = qvec*costhq

C!	
c	write(3,339) qNRS
c339	format(' NRS q-vector: ', 4f12.5)
C!
	if (Pperp .gt. 0.0) then
	  do I = 1 , 3
	  qrotNRS(I) = 0.0
	    do J = 1 , 3
	      qrotNRS(I) = qrotNRS(I) + RotMat(I,J)*qNRS(J)
	    enddo
	  enddo
	else
	  do I = 1 , 3
	    qrotNRS(I) = qNRS(I)
	  enddo
	endif

	Nu = LorMat(4,4)*NuNRS
	do I = 1 , 3
	  Nu = Nu - LorMat(4,I)*qrotNRS(I)
	enddo
c
c	Now we want to put the lab q-vector into the array ppart(1,1-4)
c
	qrotNRS(4) = NuNRS
	do I = 1 , 3
	  LorMat(4,I) = -LorMat(4,I)
	  LorMat(I,4) = -LorMat(I,4)
	enddo
	do I = 1 , 4
	  ppart(1,I) = 0.0
	  do J = 1 , 4
	    ppart(1,I) = ppart(1,I) + LorMat(I,J)*qrotNRS(J)
	  enddo
	enddo
C!	
c	write(3,338) (ppart(1,I), I = 1, 4)
c338	format(' Final q-vector: ', 4f12.5)
C!
	return
	end

	
