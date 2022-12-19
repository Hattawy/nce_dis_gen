	program Main_SimpleRad
C	
C	This program is used to generate events
C	for electrons scattering off nucleons inside
C	deuterium, INCLUDING radiative effects (taken into account only through 
C	"equivalent external radiator"
C	Based on MainDAsym.
C	The spectator momentum can be chosen.
C	Multiply by Luminosity*Acceptance*crosssection per event
C       to get total counts in each bin
C
C	Author: S.E. Kuhn	7-Feb-1994
C	Revised 13-May-1994 to allow both AO and Whitlow parametrization
C	Revised 26-Aug-2000 and 15-Feb-2001 to use up-to-date models
C	This version for inclusive asymmetries
C
C	Revised 26-Nov-2001 to generate correctly distributed events
C	Revised 23-Feb-2003 to allow generation of inclusive d events
C	Revised 11-Apr-2003 to allow light cone WF weighing
C
C	Major changes 17-7-2003 to include radiation, use elastic only.
C	Version "RadTest" to use without Keith's files
C	Use with Irad=2 in cross_Radel to get struck nucleon 4-momentum as well

C	Updated to new "newSFs" and Argonne d WF 2012

C	Re-jiggered November 2022 to replace internal radiation with eq. ext. rad.

        implicit none
	include 'instruct.inc'
	
	integer NEvents ! Limit to 100 Mio
	real*4 SigScale ! maximum reduced cross section - for singles
        real*4 plow     ! Lower cutoff for the nucleon momentum to be generated
	real*4 theplo, thephi, coslo, coshi ! range for theta (P_S rel. to beam)
        real*4 P        ! Total momentum of generated nucleon
        real*4 Px, Py, Pz       ! Cartesian Coordinates of P
	real*4 costh, sinth, phip ! Direction of P
        real*4 Weight, SigWeight ! The a priori likelihood for given limits
	real*4 Polx, Poly, Polz ! for new WF
Co	real*4 PolN
	integer mz, PolI
	real*4 mzset(3), PolPz, PolPzz, Poltest
C
	integer I, Isub, J,K,L ! separate into two loops to avoid overflow
	integer Seed
	real*4 ran1
C
C        real*4  Momentum(51), Probability(51), F11zero(51), F11piha(51),
C     *		PolConst(51) ! Arrays to hold deuteron WF information
     	real*4 NewWF(170,7) ! p [GeV], n(p), ln(int >p), P_D, TT1, TT2, TT3
        real*4 Wval, Q2, Ebeam, Hel, Sig, Q2min, Q2max, Ezero, AnglMin
        real*4 Nu, Nu0, NuMin, NuMax, DelNu, NuRat, Wlab, DelWlab, Phi, EprMin
        real*4 TanSqTh, SigP, SigN, Q2labmin
	real*4 EStar,MStar,NuNRS, Eprime, Efinal, polloss
	real*8 TotSig(150), DelSig(150) ! in bins of NU
	real*8 pvec1,pvec2,pvec3,energy
        real*4 MassN, Pi, AlPiH ! = alpha^2*4pi*hbar^2 in mbarn
	real*4 TorCur

	integer JLund, N_particles, I_Elec, I_Prot, I_Status
	real*4 ppart(4,5), cothk ! data for each particle, including struck nucleon and radiated photon
	
	real*4 elastail(131,39,240) ! param. in (E, Q^2, W)

	character*1 Input, TargType
	logical	FreeN, Prot, Newt, Lcone

        common /Constants/ Pi, MassN, AlPiH
c       common /StrFctn/ W, W1, W2, G1, G2
co	common /Datatable/ Momentum, Probability, F11zero, F11piha,
co     *          PolConst
        common /Argonne/ newWF
	common /target/ TargType
	common /randomizer/ Seed
	common /Particles/ ppart
	common /ITorus/ TorCur
	common /RadElas/ elastail

C
	SigScale = 400.0 ! Maximum value returned by Movingcross (el.)

	MassN = 0.939
	AlPiH = 2.6057E-4 ! millibarn
	Pi = 3.1415926536
Csk03	Seed = -68355
	NMC = .TRUE.
	BODEK1 = .FALSE.
	NORES = .FALSE.
	ERROR = .FALSE.
	IPOL = 1
	IPOLRES = 1
	IA1 = 4
	AsymChoice = 11
	SFChoice = 20

	do J = 1, 150
	  TotSig(J) = 0.0
	enddo

co        open (unit=31, file='deuteronwf.dat', status='OLD')
co        do J = 1, 51
co          read(31,*) Momentum(J),Probability(J), F11zero(J),
co     *		     F11piha(J), PolConst(J)
co        enddo
	open (unit=31, file='ArgonneWF.dat', status='OLD')
	read(31,*) Input
	do I = 1, 170
	  read(31,*) (NewWF(I,J), J=1,7)
	enddo
        close (unit=31)
	
	open (unit=5, file='deutin_clas12_rad.dat', status= 'OLD')
ccc	open (unit=3, file='/Volumes/nas0/user/kuhn/taggedn/4gev/deut4_radel.out', status='NEW')
	open (unit=3, file='deutout_Bonus12_eqrad.dat', status='NEW')
csk	IF particles should be written out
	open (unit=4, file='BONuS12radel_events.txt', status='NEW')
ccc	open (unit=4, file='deut4_radel_events.out', status='NEW')

csk
	write(3,*) ' Enter random number seed (+/- < 2^16):'
	read(5,*) Seed
	write(3,*) ' ', Seed
	
c	write(6,*) ' Enter STRUCK nucleon! (n/[p]):'
	write(3,*) ' Enter STRUCK nucleon! (n/[p]):'
	read (5,100) Input
	write(3,*) ' ', Input
	Prot = ((Input .ne. 'n') .and. (Input .ne. 'N'))
	Newt = ((Input .ne. 'p') .and. (Input .ne. 'P'))
100	format(A1)
csk	Newt = (.not. Prot)
csk	call TailRead (Prot) ! Read in files for radiative elastic tail

	write(3,*) ' Do you want to use light cone wave function? y[N]'
	read(5,*) Input
	write(3,*) ' ', Input
	Lcone = ((Input .eq. 'y') .or. (Input .eq. 'Y'))

c	write(6,*) ' Input Q^2 in GeV/c^2 (real number) - lower limit: '
	write(3,*) ' Input Q^2 in GeV/c^2 (real number) - lower limit: '
	read(5,*) Q2min
	write(3,*) ' ', Q2min
	write(3,*) ' Input Q^2 in GeV/c^2 (real number) - lower LAB limit (after radiation): '
	read(5,*) Q2labmin
	write(3,*) ' ', Q2labmin
c	write(6,*) ' Input Q^2 in GeV/c^2 (real number) - upper limit: '
	write(3,*) ' Input Q^2 in GeV/c^2 (real number) - upper limit: '
	read(5,*) Q2max
	write(3,*) ' ', Q2max
	
csk	if (Prot) SigScale = Sigscale*2.0
	SigScale = SigScale/(1.0+Q2min/0.71)**4
	if (SigScale .lt. 1.0) SigScale = 1.0
	if (Lcone) SigScale = SigScale*2.0 ! LC correction is 0-2
	if (Prot.and.Newt) SigScale = SigScale*2.0

c	write(6,*) ' Input Beam Energy: '
	write(3,*) ' Input Beam Energy: '
	read(5,*) Ezero
	write(3,*) ' ', Ezero	
	
911	write(3,*) ' Input nu in GeV (real number) - upper limit: '
c	write(6,*) ' Input nu in GeV (real number) - upper limit: '
	read(5,*) NuMax
csk	if (NuMax .gt. Ebeam) goto 911
	write(3,*) ' ', NuMax
	
912	write(3,*) ' Input nu in GeV (real number) - lower limit: '
c	write(6,*) ' Input nu in GeV (real number) - lower limit: '
	read(5,*) NuMin
	if (NuMin .ge. NuMax) goto 912
	DelNu = NuMax-NuMin
	NuRat = alog(NuMax/NuMin)
	write(3,*) ' ', NuMin

c	write(6,*) ' Non-resonant part only? y/[N]'
	write(3,*) ' Non-resonant part only? y/[N]'
	read (5,100) Input
	write(3,*) ' ', Input
	NORES = ((Input .eq. 'y') .or. (Input .eq. 'Y'))
	
	if (Newt .and. Prot) then
	  write(3,*) ' Inclusive scattering on deuterium.'
c	  write(6,*) ' Inclusive scattering on deuterium.'
	endif
	  
913	write(3,*) ' Lower limit for spectator momentum: '
c	write(6,*) ' Lower limit for spectator momentum: '
	read(5,*) plow
	if (plow .ge. 0.7) then
	  write(6,*) ' Must be below 0.7 GeV!'
	  goto 913
	else if (plow .lt. 0.0) then
	  FreeN = .TRUE.	! Assume free N cross section
c	  write(6,*) ' Assume free nucleon! '
	  write(3,*) ' Assume free nucleon! '
	else
	  FreeN = .FALSE.
	endif
	write(3,*) ' ', plow

13	write(3,*) ' Most backward spectator angle (Lab, wrt beam): '
c	write(6,*) ' Most backward spectator angle (Lab, wrt beam): '
	read(5,*) thephi
	write(3,*) ' ', thephi
	write(3,*) ' Most forward spectator angle: '
c	write(6,*) ' Most forward spectator angle: '
	read(5,*) theplo
	if (thephi .le. theplo) goto 13
	write(3,*) ' ', theplo
	coslo = cos(theplo*pi/180.0)
	coshi = cos(thephi*pi/180.0)

1214	write(3,*)
     *	 ' Enter Pz (-1...1) and Pzz (-2...1) of the deuteron NOW: '
c	write(6,*)
c     *	 ' Enter Pz (-1...1) and Pzz (-2...1) of the deuteron NOW: '
	read(5,*) PolPz
	read(5,*) PolPzz
	PolPzz = (PolPzz + 2.0) / 3.0
	mzset(1) = (PolPz + PolPzz) / 2.0
	mzset(3) = (PolPzz - PolPz) / 2.0
	mzset(2) = 1.0 - mzset(1) - mzset(3)
	write(3,*) 'N+ , N0, N-: ',mzset
c	write(6,*) 'N+ , N0, N-: ',mzset
	if ((mzset(1) .lt. 0.0) .or. (mzset(2) .lt. 0.0) .or.
     2		(mzset(3) .lt. 0.0)) goto 1214

	mz = 0
	if (FreeN) then
	  Weight = 1.0
	else
co	  call Deuteron (plow,mz,coslo,coshi, P,Px,Py,Pz,PolN,Weight)
	  call deuteron_new (plow,mz,coslo,coshi, P,Px,Py,Pz,Polx,Poly,Polz,Weight)
	endif
	write(3,*) ' Weight for selection of d WF: ', Weight
	SigWeight = AlPiH*(1/Q2min - 1/Q2max)*NuRat
	write(3,*) ' Weight for cross section: ', SigWeight
	write(3,*) ' Maximum Scale for reduced Sigma: ', SigScale
	
	write(3,*) ' How many simulated events?'
c	write(6,*) ' How many simulated events?'
	read(5,*) NEvents
	write(3,*) ' ', NEvents
	SigWeight = Weight*SigWeight/Nevents/2000.*SigScale
C	/2 for hel+/-, SigScale for max. reduced cross section
	write(3,*) ' Cross section for each event: [millibarn] ',  SigWeight
	
	write(3,*) ' Torus current? (up to +/- 3375)'
c	write(6,*) ' Torus current? (up to +/- 3375)'
	read(5,*) TorCur
	write(3,*) ' ', TorCur, ' Amps'
	write(3,*) ' Minimum Lab scattering angle for e-? [deg]'
c	write(6,*) ' Minimum Lab scattering angle for e-? [deg]'
	read(5,*) AnglMin
	write(3,*) ' ', AnglMin, ' degrees'
	AnglMin = 4.0*((sin(AnglMin*Pi/360.))**2)
	write(3,*) ' Minimum Lab Eprime? [GeV]'
c	write(6,*) ' Minimum Lab Eprime? [GeV]'
	read(5,*) EprMin
	write(3,*) ' ', EprMin, ' GeV'
c	write(6,*) ' Starting Random Nr.:',ran1(Seed)
	write(3,*) ' Starting Random Nr.:',ran1(Seed)
	close (unit=5)
		
        do I = 1 , NEvents
	  if((I/100000)*100000.eq.I) write(3,*) I
	do Isub = 1, 1000 ! nested loop to avoid overflows

	  Q2 = Q2min/(1.0 - ran1(Seed)*(1.0 - Q2min/Q2max))
	 
	  call radbefB12(Seed,Q2,Ezero,Ebeam,polloss)
cssk	  polloss = 0.0 ! cssk
cssk	  Ebeam = Ezero ! cssk
	  if (Ebeam .le. EprMin) goto 113
	  if(Q2.lt.(Ebeam*EprMin*AnglMin)) goto 113
	  
	  Poltest = ran1(Seed)
	  if (Poltest .le. mzset(1)) then
	    mz = +1
	  else if (Poltest .le. (mzset(1) + mzset(2))) then
	    mz = 0
	  else
	    mz = -1
	  endif
	  if (FreeN) then
	    P = 0.0
	    Px = 0.0
	    Py = 0.0
	    Pz = 0.0
	    Polx = 0.0
	    Poly = 0.0
	    Polz = mz
	    Weight = 1.0
	  else
co14          call Deuteron (plow,mz,coslo,coshi, P,Px,Py,Pz,PolN,Weight)
	  call deuteron_new (plow,mz,coslo,coshi, P,Px,Py,Pz,Polx,Poly,Polz,Weight)
            if (P .ge. 0.7) goto 113
	    
	   Px = -Px
	   Py = -Py
	   Pz = -Pz ! I believe that deuteron_new produces the spectator mom. -> need struck N
C          If P>0.7, MStar becomes imaginary
	  endif
	  Polx = Polx*(1. - polloss)
	  Poly = Poly*(1. - polloss)
	  Polz = Polz*(1. - polloss)

          EStar = 1.876 - sqrt(MassN*MassN + P*P)
	  if (FreeN) Estar = MassN
          MStar = sqrt(EStar*EStar - P*P)
csk        NuNRS = DelNu*ran1(Seed) + NuMin
	  NUNRS = NuMin*exp(NuRat*ran1(Seed))
	  Phi = 2*pi*ran1(Seed) ! angle of scattered electron
	  cothk=2.*ran1(Seed) - 1.0 ! angle of radiated photon - probably not needed

          if ((MStar*(MStar+2.0*NuNRS)) .gt. (0.87 + Q2)) then ! W big enuff
C
C	First for positive helicity
C
co            Hel = PolN
	    Hel = 1.0
	    Sig = 0.0
	    if (Prot) then
	      TargType = 'P'
               call MovingCross_new (P, Px, Py, Pz, MStar, EStar,
     *            Q2, Ebeam, Hel, Polx, Poly, Polz, NuNRS, Phi, SigP, Nu0)
	        Sig = SigP
	    endif
	    if (Newt) then
	      TargType = 'N'
co              call MoCrossRadel (P, Px, Py, Pz, MStar, EStar,
co     *            Q2, Ebeam, Hel, NuNRS, Phi, SigN, Nu0,cothk)
               call MovingCross_new (P, Px, Py, Pz, MStar, EStar,
     *            Q2, Ebeam, Hel, Polx, Poly, Polz, NuNRS, Phi, SigN, Nu0)
     	      Sig = Sig + SigN
	    endif
	    if ((Nu0 .le. 0.0) .or. (Nu0 .ge. Ebeam-EprMin)) goto 112
	    Sig = Sig/SigScale*NuNRS
	    Eprime = Ebeam - Nu0	! Scattered Electron
	    call radaftB12(Seed,Q2,Eprime,Efinal,polloss) ! cludge for after-scattering radiative loss
cssk	    Efinal = Eprime !cssk
	    
	    if (Sig .lt. 1.0e-11) goto 112
	    if (Efinal .lt. EprMin) goto 112
	    
csk	    Struck nucleon 4-momentum - needed in writeout
	    ppart(3,1) = Px + ppart(1,1)
	    ppart(3,2) = Py + ppart(1,2)
	    ppart(3,3) = Pz + ppart(1,3)
	    ppart(3,4) = Estar + Nu0
	    
	    ppart(1,4) = Efinal
csk		NEW! 
	    ppart(4,4) = Eprime - Efinal ! radiated photon
	    ppart(1,1) = -ppart(1,1)*Efinal/Eprime
	    ppart(4,1) = ppart(1,1)*ppart(4,4)/Efinal
	    ppart(1,2) = -ppart(1,2)*Efinal/Eprime
	    ppart(4,2) = ppart(1,2)*ppart(4,4)/Efinal
	    ppart(1,3) = (Ebeam - ppart(1,3))*Efinal/Eprime ! based on q-vector
	    ppart(4,3) = ppart(1,3)*ppart(4,4)/Efinal
	    ppart(4,5) = 0.0
	    
	    if((ppart(1,3)/Efinal) .gt. (1.-AnglMin/2.)) goto 112
	    if ((2*Ezero*(Efinal-ppart(1,3))) .lt. Q2labmin) goto 112

	    ppart(1,5) = 0.000511
	    ppart(2,1) = -Px	! Recoil nucleon
	    ppart(2,2) = -Py
	    ppart(2,3) = -Pz
	    ppart(2,4) = 1.876 - EStar
	    ppart(2,5) = MassN
	    
	    ppart(3,5) = MassN ! Not used here

c	    ppart(1,4) = Ebeam - Nu0	! Scattered Electron
c	    ppart(1,1) = -ppart(1,1)
c	    ppart(1,2) = -ppart(1,2)
c	    ppart(1,3) = Ebeam - ppart(1,3) ! based on q-vector
c	    ppart(1,5) = 0.000511
c	    ppart(2,1) = -Px	! Recoil nucleon
c	    ppart(2,2) = -Py
c	    ppart(2,3) = -Pz
c	    ppart(2,4) = 1.876 - EStar
c	    ppart(2,5) = MassN

	    call writeout(1, SigWeight, Sig, Ezero, Prot, Newt, Lcone, 1)
	    
111	    continue
C
C	Now for negative helicity
C
co            Hel = -PolN
	    Hel = -1.0
	    Sig = 0.0
	    if (Prot) then
	      TargType = 'P'
co              call MoCrossRadel (P, Px, Py, Pz, MStar, EStar,
co     *            Q2, Ebeam, Hel, NuNRS, Phi, SigP, Nu0,cothk)
               call MovingCross_new (P, Px, Py, Pz, MStar, EStar,
     *            Q2, Ebeam, Hel, Polx, Poly, Polz, NuNRS, Phi, SigP, Nu0)
     	      Sig = SigP
	    endif
	    if (Newt) then
	      TargType = 'N'
co              call MoCrossRadel (P, Px, Py, Pz, MStar, EStar,
co     *            Q2, Ebeam, Hel, NuNRS, Phi, SigN, Nu0,cothk)
               call MovingCross_new (P, Px, Py, Pz, MStar, EStar,
     *            Q2, Ebeam, Hel, Polx, Poly, Polz, NuNRS, Phi, SigN, Nu0)
     	      Sig = Sig + SigN
	    endif

csk	    if ((Nu0 .le. 0.0) .or. (Nu0 .ge. Ebeam)) goto 112
	    Sig = Sig/SigScale*NuNRS
c	    if (Sig .lt. 1.0e-11) goto 112
	    
csk	    Struck nucleon 4-momentum - needed in writeout
	    ppart(3,1) = Px + ppart(1,1)
	    ppart(3,2) = Py + ppart(1,2)
	    ppart(3,3) = Pz + ppart(1,3)
	    ppart(3,4) = Estar + Nu0
	    
	    ppart(1,4) = Efinal
csk		NEW! 
	    ppart(4,4) = Eprime - Efinal ! radiated photon
	    ppart(1,1) = -ppart(1,1)*Efinal/Eprime
	    ppart(4,1) = ppart(1,1)*ppart(4,4)/Efinal
	    ppart(1,2) = -ppart(1,2)*Efinal/Eprime
	    ppart(4,2) = ppart(1,2)*ppart(4,4)/Efinal
	    ppart(1,3) = (Ebeam - ppart(1,3))*Efinal/Eprime ! based on q-vector
	    ppart(4,3) = ppart(1,3)*ppart(4,4)/Efinal
	    ppart(4,5) = 0.0
	    ppart(1,5) = 0.000511
	    ppart(2,1) = -Px	! Recoil nucleon
	    ppart(2,2) = -Py
	    ppart(2,3) = -Pz
	    ppart(2,4) = 1.876 - EStar
	    ppart(2,5) = MassN
	    
	    ppart(3,5) = MassN

c	    ppart(1,4) = Ebeam - Nu0	! Scattered Electron
c	    ppart(1,1) = -ppart(1,1)
c	    ppart(1,2) = -ppart(1,2)
c	    ppart(1,3) = Ebeam - ppart(1,3) ! based on q-vector
c	    ppart(1,5) = 0.000511
c	    ppart(2,1) = -Px	! Recoil Proton
c	    ppart(2,2) = -Py
c	    ppart(2,3) = -Pz
c	    ppart(2,4) = 1.876 - EStar
c	    ppart(2,5) = MassN
	      
	    call writeout(1, SigWeight, Sig, Ezero, Prot, Newt, Lcone, 2)
	    
112	    continue
          endif ! W big enuff
113	  continue
	enddo ! inner of two nested loops to avoid overflow
	enddo
        close(unit=4)
csk	 if (.not. Prot) then
	call writeout(2, SigWeight, Sig, Ezero, Prot, Newt, Lcone, 3)
csk	 endif
	close(unit=3)
	stop
	end
