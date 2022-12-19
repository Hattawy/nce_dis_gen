
	subroutine writeout(I, SigWeight, Sig, E, Prot, Newt, Lcone, Ihel)
	
C	Subroutine to fill arrays with information about a MC run and print out events in lund format

	implicit none
	real*4 ppart(4,5), pmiss(4)	! data for each particle: e, true spec, struck N, reconstr. spec
	real*4 Estrue, Pstrue, psqtrue, cosPsqtrue, alphastrue,
     *		Mstartrue, Estartrue, Wstartrue ! Variables derived from real spectator momentum ppart(2,i)
        real*4 qtrue(3), qvaltrue ! We have to infer the real three-momentum transfer
        real*4 lcweight, Hel, SigWeight, Sig, SigP, SigCheck
	real*8 Weight
	integer I, J, K, L, MM, Ihel ! Don't use I, Ihel - already taken!!!
	integer IQ, Ip, Icos, Iphi, ISector, id ! particle ID
	integer Seed
	real*4 Pi /3.1415926536/
	real*4 E, Torus
	real*4 pacc, thetacc, phiacc, toracc, dphiacc, weightacc, eldelphi ! for FASTMC
		
	real*4 qvect(4), qval, nu, Q2, thete, phie, xBjLab, xilab, Wlab ! Electron-dependent variables
	real*4 PmissT, pmissq, cosPmissq, alphamiss, pmiss_perp, pmiss_long, emiss, mmiss, 
     *		Mstar, Estar, Wstar, xBjnrs, xinrs !spectator-dependent variables used for binning
	real*4 Pft(4), Pftot, Pfq, cosPf, thetp, phip ! struck nucleon if it is intact in final state
	real*4 MassN, MassD
	real*4 cx, cy, sinthemin, phisector
	real*4 expon, delphi
	real*8 nudis(100), nudis1cut(100), nudis2cut(100), nudis3cut(100)
	real*8 Q2dis(150,2), Q2dis1cut(150,2), Q2dis2cut(150,2), Q2dis3cut(150,2)
	real*8 Psdis(100), Psdis1cut(100), Psdis2cut(100), Psdis3cut(100)
	real*8 thpqdis(90), thdis1cut(90), thdis2cut(90), thdis3cut(90)
	real*8 Wstardis(450,2), Wsdis1cut(450,2), Wsdis2cut(450,2),
     *	 Wsdis3cut(450,2)
	real*8 Wlabdis(450,2), Wldis1cut(450,2), Wldis2cut(450,2),
     *	 Wldis3cut(450,2)
	real*8 nuQ2dis(100,30), xialpdis(100,20), ximstdis(100,20)
	real*8 WQ2Sumdis(450,45), WQ2Diffdis(450,45) ! W=0-4.5, 45 logarithmic Q2 bins
	real*8 xQ2dis(20,15), xQ2dis2(20,15) ! 0.05 steps in x, 1 GeV/c2 steps in Q2
	real*8 mikeDis(40,5,3,2), mikeDis2(40,5,9,2)! IQ;  Ip, Icos, Ihel
	real*4 ran1
     	real*4 NewWF(170,7) ! p [GeV], n(p), ln(int >p), P_D, TT1, TT2, TT3
	real*4 AA /41./, BB /0.26/, CC /0.3/, DD(2) /9.,8./, EE /16.72/, FF /0.06/, apr, bpr
	real*4 Philimit /10.0/, ThetaHigh /45.0/, themin /5.0/
	integer imombin
	
csk	Variables needed for lund file output including accidentals
	real*8 bgndweight(34,20) ! 34 ps bins, 20 costheta
	real*8 vx/0.0/, vy/0.0/, vz, Mel, mph/0.0/, wiwi, zero/0.0/
	real*8 pxpbg, pypbg, pzpbg, Epbg, phibg, costhbg,ppbg,eventweight/1.0/
	integer Npart,npfinal,Atarg/2/, nbgnd/0/
	integer Ipp/2212/,In/2112/,Iel/11/,Iph/22/,null/0/,one/1/,irow
	logical realp/.TRUE./
	
	logical first /.TRUE./, Cut1, Cut2, Cut3
	logical	Prot, Newt, Lcone, iselast
	
	common /randomizer/ Seed
	common /ITorus/ Torus
	common /Particles/ ppart
        common /Argonne/ newWF
	common /Safer/ nudis, nudis1cut, nudis2cut, nudis3cut,
     *		Q2dis, Q2dis1cut, Q2dis2cut, Q2dis3cut,
     *		Psdis, Psdis1cut, Psdis2cut, Psdis3cut,
     *		thpqdis, thdis1cut, thdis2cut, thdis3cut,
     *		Wstardis, Wsdis1cut, Wsdis2cut, Wsdis3cut,
     *		Wlabdis, Wldis1cut, Wldis2cut, Wldis3cut,
     *		nuQ2dis, xialpdis, ximstdis, WQ2Sumdis, WQ2Diffdis,xQ2dis,
     *		mikeDis, mikeDis2, bgndweight

	MassN = 0.939
	MassD = 1.876 ! Deuteron Hardcode !
	Mel = 0.000511
	if (first) then
	  print *, ' First call of writeout.f'
	  first = .FALSE.
	  do IQ = 1, 40
	    do Ip = 1, 5
	      do J = 1,2
csk	      do Iphi = 1,4
	        do Icos = 1, 3
		  mikeDis(IQ,Ip,Icos,J) = 0.0D0
		enddo
		do Icos = 1,9
		  mikeDis2(IQ,Ip,Icos,J) = 0.0D0
		enddo
	      enddo
	    enddo
	  enddo

	  do J = 1, 100
	    nudis(J) = 0.0
	    nudis1cut(J) = 0.0
	    nudis2cut(J) = 0.0
	    nudis3cut(J) = 0.0
	    Psdis(J) = 0.0
	    Psdis1cut(J) = 0.0
	    Psdis2cut(J) = 0.0
	    Psdis3cut(J) = 0.0
	  enddo
	  do J = 1, 150
	    do Ip = 1, 2
	    Q2dis(J,Ip) = 0.0
	    Q2dis1cut(J,Ip) = 0.0
	    Q2dis2cut(J,Ip) = 0.0
	    Q2dis3cut(J,Ip) = 0.0
	    enddo
	  enddo
	  do J = 1, 90
	    thpqdis(J) = 0.0
	    thdis1cut(J) = 0.0
	    thdis2cut(J) = 0.0
	    thdis3cut(J) = 0.0
	  enddo
	  do J = 1, 450
	  do L = 1, 2
	    Wstardis(J,L) = 0.0
	    Wsdis1cut(J,L) = 0.0
	    Wsdis2cut(J,L) = 0.0
	    Wsdis3cut(J,L) = 0.0
	    Wlabdis(J,L) = 0.0
	    Wldis1cut(J,L) = 0.0
	    Wldis2cut(J,L) = 0.0
	    Wldis3cut(J,L) = 0.0
	  enddo
	  enddo
	  do J = 1, 100
	    do L = 1, 20
	      xialpdis(J,L) = 0.0
	      ximstdis(J,L) = 0.0
	    enddo
	    do L = 1, 30
	      nuQ2dis(J,L) = 0.0
	    enddo
	  enddo
	  do J = 1, 450
	    do L = 1, 45
	      WQ2Sumdis(J,L) = 0.0
	      WQ2Diffdis(J,L) = 0.0
	    enddo
	  enddo
          do J = 1, 20
            do L = 1, 15
              xQ2dis(J,L) = 0.0
              xQ2dis2(J,L) = 0.0
            enddo
          enddo

c	Integrated cross section file for background protons in 20 costheta bins
c	and 34 ps pins from 0.06-0.07 to 0.39-0.4 GeV/c. Used to write background
csk	  do L = 1,34
c	    do J = 1,20
c	      read(23,*) bgndweight(L,J)
c	    enddo
c	  enddo
c	
	  Npart=nbgnd+1 ! + electron
	  if(realp)Npart = Npart+1 ! +Spectator
c	  print *, ' Made it through initialization', Npart
	endif ! "first"

	if (I .eq. 1) then
	if (ppart(1,3) .le. 0.001) return
	  qvect(1) = -ppart(1,1)
	  qvect(2) = -ppart(1,2)
	  qvect(3) = E - ppart(1,3)
	  qvect(4) = E - ppart(1,4)
	  cx = ppart(1,1)/ppart(1,4)
	  cy = ppart(1,2)/ppart(1,4)
	  thete = 180./pi*acos(ppart(1,3)/ppart(1,4))
5827	format(' ',8f10.6)
	  
	  phie = 180./Pi*atan2(cy,cx)
C
	  phisector = phie
C
	  if (phisector .lt. -30.0) phisector = phisector + 360.0
c	  phie = phisector + 30.0 ! Only for debugging purposes
	  
	  if (phisector .gt. 30.0) then
	    Isector = (phisector + 30.0)/60.0
	    phisector = phisector - 60.0*Isector
	  endif
	  
9991	  format(' ',4f10.4)
	  nu = qvect(4)
	  qval = sqrt(qvect(1)**2+qvect(2)**2+qvect(3)**2)
	  Q2 = qval*qval - nu*nu
CCCC	  write(3,*) ' RECON Q2:', Q2
	  xBjLab = Q2/2/MassN/nu
	  Wlab = sqrt(MassN*MassN + 2.0*MassN*nu - Q2)
	  xilab = (qval - nu)/MassN

C	True spectator values needed for LC correction ! Deuteron Hardcode !
	  Estrue = ppart(2,4)
	  Pstrue = sqrt(ppart(2,1)**2+ppart(2,2)**2+ppart(2,3)**2)
	  pmissq = 0.0 ! pmissq is the reconstructed product incl. radiative effects
	  do J = 1, 3
	    pmissq = pmissq + qvect(J)*ppart(2,J)
	  enddo
	  
	  qvaltrue = 0.0 ! Here we calculate true lc variables
	  psqtrue = 0.0
	  do J = 1 , 3
	    qtrue(J) = ppart(3,J) + ppart(2,J)
	    psqtrue = psqtrue + qtrue(J)*ppart(2,J)
	    qvaltrue = qvaltrue + qtrue(J)*qtrue(J)
	  enddo
	    qvaltrue = sqrt(qvaltrue)
     	  if ((Pstrue*qvaltrue) .gt. 0.0) then
	    alphastrue = (Estrue - psqtrue/qvaltrue)/MassN
	    cosPsqtrue = psqtrue/qvaltrue/Pstrue
	  else
	    cosPsqtrue = 0.0
	    alphastrue = 1.0
	  endif
	  if (alphastrue .ge. 2.0) alphastrue = 1.9999
	  if (alphastrue .le. 0.0) alphastrue = 0.0001
	  if  (abs(cosPsqtrue).ge. 1.0) cosPsqtrue = cosPsqtrue/abs(cosPsqtrue)*0.9999
	  Estartrue = MassD - Estrue
	  Mstartrue = sqrt(Estartrue*Estartrue-Pstrue*Pstrue)

	  SigCheck = Sig
	  if (Lcone) then
	    call LCsupp_new (Pstrue,alphastrue, cosPsqtrue, lcweight)
	    SigCheck = Sig*lcweight
	    SigP = SigP*lcweight
	  endif
c	  print *, ' Made it through Lcone'
	  
CCC	The following is if we can USE the spectator momentum in ppart(2,i)	  
	  
	  PmissT = Pstrue
	  pmissq = (qvect(1)*ppart(2,1)+qvect(2)*ppart(2,2)+
     *		qvect(3)*ppart(2,3))
	  emiss = Estrue
	  mmiss = MassN
	  Estar = MassD - emiss
	  pmiss_long = ppart(2,3)
	  pmiss_perp = sqrt(ppart(2,1)**2 + ppart(2,2)**2)
	  
CCC	Alternatively, here we assume the struck nucleon was detected in ppart(3,i)
	  
c	  PmissT = sqrt((ppart(3,1)-qvect(1))**2+(ppart(3,2)-qvect(2))**2+
c     & (ppart(3,3)-qvect(3))**2)
c	  pmissq = (qvect(1)*ppart(3,1)+qvect(2)*ppart(3,2)+
c     *		qvect(3)*ppart(3,3))-qval*qval
c	  emiss = MassD + nu - ppart(3,4) ! Deuteron Hardcode !
c	  mmiss = sqrt(emiss*emiss - PmissT*PmissT)
c	  Estar = MassD - emiss ! Deuteron Hardcode !
c	  pmiss_long = ppart(3,3) - qvect(3)
c	  pmiss_perp = sqrt(PmissT*PmissT - pmiss_long*pmiss_long)
	  
CCC
	  Mstar = sqrt(Estar*Estar-PmissT*PmissT) ! Deuteron Hardcode !
	  Wstar = sqrt(Mstar*Mstar + 2.0*(Estar*nu+pmissq)-Q2) ! Deuteron Hardcode !
c          Wstar = ppart(3,5)
	  iselast = ((Wstar .lt. 1.00).and.(Wstar .gt. 0.9))

ccc	  if (iselast) then
ccc	    do J = 1, 3
ccc	      Pft(J) = qvect(J) - ppart(2,J)
ccc	    enddo
ccc	    Pftot = sqrt(Pft(1)**2 + Pft(2)**2 + Pft(3)**2)
ccc	    thetp = 180./pi*acos(Pft(3)/Pftot)
ccc	    phip = 180./Pi*atan2(Pft(2),Pft(1))
ccc	    if (phip .lt. -30.0) phip = phip + 360.0
ccc	    if (phip .gt. 30.0) then
ccc	      Isector = (phip + 30.0)/60.0
ccc	      phip = phip - 60.0*Isector
ccc	    endif
	      
ccc	    Pfq = qval**2 - pmissq
ccc     	    if ((Pftot*qval) .gt. 0.0) then
ccc	      cosPf = Pfq/qval/Pftot
ccc	    else
ccc	      cosPf = 0.0
ccc	    endif
ccc	  endif	  
	  
          if ((PmissT*qval) .gt. 0.0) then
	    cosPmissq = pmissq/qval/PmissT
	    alphamiss = (emiss - pmissq/qval)/MassN ! Deuteron Hardcode !
	  else
	    cosPmissq = 0.0
	    alphamiss = 1.0
	  endif
	  if (alphamiss .ge. 2.0) alphamiss = 1.9999
	  if (alphamiss .le. 0.0) alphamiss = 0.0001	  
CCC
	  if (alphamiss .lt. 2.0) then ! Deuteron Hardcode !
	    xinrs = xilab/(2.0 - alphamiss)
	  else
	    xinrs = 1.0
	  endif
	  xBjnrs = Q2/2/(Estar*nu+pmissq) ! Deuteron Hardcode !
	  
c	print *, ' All but cuts are calculated'
	  Cut1 = ((ppart(1,4).ge.0.1*E).and.(ppart(1,4) .lt. E)) ! used to be 0.5!
c	  Cut1 = Cut1 .and.((thete.gt.themin).and.(thete.le.ThetaHigh))
c	  if (Cut1) then
c	    id = 11
c	    pacc = ppart(1,4)
c	    thetacc = thete
c	    phiacc = phisector
c	    toracc = Torus
c	    call clas_at12g(id, pacc, thetacc, phiacc, toracc, dphiacc, weightacc)
c	    Cut1 = (dphiacc .gt. 0.0)
c	    eldelphi = dphiacc
c	  endif
c
c	  Cut2 = (Cut1 .and. weightacc .ge. 0.9)
	  Cut2 = Cut1 .and.((thete.gt.themin).and.(thete.le.ThetaHigh))
	  Cut3 = Cut2
c	  if(iselast.and.Cut2 .and. SigP .gt. 0.0) then
c	    id = 2212
c	    pacc = sqrt(ppart(3,1)**2 + ppart(3,2)**2 + ppart(3,3)**2)
c	    thetacc = 180./pi*acos(ppart(3,3)/pacc)
c	    phiacc = 180./Pi*atan2(ppart(3,2),ppart(3,1))
c	    if (phiacc .lt. -30.0) phiacc = phiacc + 360.0
c	    if (phiacc .gt. 30.0) then
c	      Isector = (phiacc + 30.0)/60.0
c	      phiacc = phiacc - 60.0*Isector
c	    endif
c	    call clas_at12g(id, pacc, thetacc, phiacc, toracc, dphiacc, weightacc)
c	    Cut3 = ((weightacc .gt. 0.9) .and. (PmissT.lt.0.3)) !Cludge for 7Li
c	  endif
c      if (Cut2) then
c	    Cut3 = ((PmissT.lt.0.5).and.(emiss.lt.0.15).and.(pmiss_long.lt.0.2))
c      endif
CSK       Here goes the regular printout of singles
                    
          Weight = SigWeight*SigCheck
c	  print *, ' here we are RIGHT before writeout'
CSK	  goto 111 ! Comment out for writeout
	  if(.not.Cut3) goto 111 
CSK writeout           
CSK	  do K = 1, 1000 ! For the rare occasion that SigCheck>1
CSK         if (SigCheck .lt. ran1(Seed)) goto 111
CSK	    if (.not. (Prot.and.Newt)) then
CSK	      write(4,1221) (ppart(1,J), J=1,4),(ppart(2,J), J=1,4),(ppart(3,J), J=1,4)
CSK	    else
CSK	      write(4,1020) (ppart(1,J), J=1,4)
CSK	    endif
CSK	    SigCheck = SigCheck - 1.0
CSK	  enddo
CSK1020	  format(' ',4f8.4)
CSK1221	  format(' ',12f8.4)
CSK	  if (SigCheck .gt. 0.0) then
CSK	    write(6,*) ' SigCheck too large after 1000 iterations', SigCheck, Sig, SigWeight, Weight
CSK	  endif
CSK end writeout 

CCC	Writeout of complete Lund format with background
	do K = 1, 1000 ! For the rare occasion that SigCheck>1
	  if (SigCheck .lt. ran1(Seed)) goto 111
c	print *, ' Begin writeout'
	  npfinal = Npart
	  if(ppart(4,4).gt.0.05)npfinal = Npart+1
	  write(4,101) npfinal,Atarg,Atarg,zero,zero,Iel,E,In,one,eventweight ! Event header
101	format(' ',I2,2I4,2f5.1,I5,f8.3,2I5,E12.3)	
	  vz = 40.0*ran1(Seed)-20.0
	  irow = 1 ! scattered electron
	  write(4,102) irow,zero,one,Iel,null,null,(ppart(1,J), J=1,4),mel,vx,vy,vz
102	format(' ',I2,f4.0,I4,3I5,8f10.4)
	  if(npfinal.gt.Npart)then
	    irow = irow+1 ! radiated photon
	    write(4,102) irow,zero,one,Iph,null,null,(ppart(4,J), J=1,4),mph,vx,vy,vz
	  endif	  
	  if(realp)then
	    irow = irow+1
	    write(4,102) irow,zero,one,Ipp,null,null,(ppart(2,J), J=1,4),MassN,vx,vy,vz 
	  endif
c	  print *, ' Wrote e and real p'
	  
	  if (nbgnd .gt. 0) then
	  do J = 1, nbgnd
	    irow = irow+1
	    wiwi = ran1(Seed)
	    ppbg = 0.06
	    do L = 1,34
	    costhbg = -0.99999
	      do MM = 1,20
	        if (bgndweight(L,MM) .ge. wiwi) then
		  ppbg = ppbg + ran1(Seed)*0.01
		  costhbg = costhbg + ran1(Seed)*0.09999
		  goto 991
		endif
		costhbg = costhbg + 0.1
	      enddo
	    ppbg = ppbg + 0.01
	    enddo
991	    continue ! found a randomly distributed spectator proton
	    pzpbg = ppbg*costhbg
	    pxpbg = ppbg*dsqrt(1.0D0-costhbg*costhbg)
	    phibg = 2.0D0*pi*ran1(Seed)
	    pypbg = pxpbg*dsin(phibg)
	    pxpbg = pxpbg*dcos(phibg)
	    Epbg = dsqrt(MassN*MassN + ppbg*ppbg)
	    vz = 40.0*ran1(Seed)-20.0
	    write(4,102) irow,zero,one,Ipp,null,null,pxpbg,pypbg,pzpbg,Epbg,MassN,vx,vy,vz 
	  enddo ! finished writing background protons
	  endif
	SigCheck = SigCHeck - 1.0
	enddo
ccc	END LUND WRITEOUT
          
111	  continue
c	print *, ' Start histogramming'
c	if (Q2 .lt. 0.13 .or. Q2 .gt. 3.18) goto 1112
c	L = 27+13.0*(alog(Q2/0.91883882)/alog(10.0)) ! EG1b Q2 bins
c	if (L.lt.1 .or. L.gt.40) goto 1112
c	Ip = PmissT*10.0+1.5
c	if (Ip .gt. 5 .and. PmissT.lt.0.5) Ip = 5
c	if (Ip .lt. 1 .or. Ip .gt. 5) goto 1112
c	Icos = 1
c	if(cosPmissq .gt. -0.35) Icos = 2
c	if(cosPmissq .gt. 0.35) Icos = 3
c	if (Icos .lt. 1 .or. Icos .gt. 3) goto 1112
cc	Iphi = 1
cc	if (Iphi .lt. 1 .or. Iphi .gt. 4) goto 1112
c	if (Cut3) then
c	    mikeDis(L,Ip,Icos,1) = mikeDis(L,Ip,Icos,1)  + Weight
c	    mikeDis(L,Ip,Icos,2) = mikeDis(L,Ip,Icos,2)  + Weight*(3.0-2.0*Ihel)
c	endif	    

1112	  K = (Nu/0.1 + 1.0)
	  if(K.gt.0 .and. K.le. 100) then
	    nudis(K) = nudis(K) + Weight
	    if(Cut1) nudis1cut(K) = nudis1cut(K) + Weight!*eldelphi
	    if(Cut2) nudis2cut(K) = nudis2cut(K) + Weight
	    if(Cut3) nudis3cut(K) = nudis3cut(K) + Weight!*SigP
	    L = (Q2/0.5 + 1.0)
	    if((L.gt.0 .and. L.le.30).and. Cut1) then
	    nuQ2dis(K,L) = nuQ2dis(K,L) + Weight
	    endif
	  endif
c	  print *, ' First Histo'
	  L = (Q2/0.1 + 1.0)
	  if(L.gt.0 .and. L.le. 150) then
c	  print *, ' second histo: L, Ihel', L, Ihel, Cut1, Cut2, Cut3
	    Q2dis(L,Ihel) = Q2dis(L,Ihel) + Weight
	    if(Cut1) Q2dis1cut(L,Ihel) = Q2dis1cut(L,Ihel) + Weight!*eldelphi
	    if(Cut2) Q2dis2cut(L,IHel) = Q2dis2cut(L,IHel) + Weight
	    if(Cut3) Q2dis3cut(L,IHel) = Q2dis3cut(L,IHel) + Weight!*SigP
	  endif
c	  print *, ' Second histo done'
	  K = (PmissT/0.01 + 1.0)
csk	  K = (Pftot/0.1 + 1.0)
	  if(K.gt.0 .and. K.le. 100) then
csk	  if(K.gt.0 .and. K.le. 100 .and. iselast) then
	    Psdis(K) = Psdis(K) + Weight
	    if(Cut1) Psdis1cut(K) = Psdis1cut(K) + Weight!*eldelphi
	    if(Cut2) Psdis2cut(K) = Psdis2cut(K) + Weight
	    if(Cut3) Psdis3cut(K) = Psdis3cut(K) + Weight!*SigP
	  endif
c	  print *, ' Third'
c	  K = ((acos(cosPmissq)*90/pi) + 1.0) ! measured theta_pq
	  K = ((acos(cosPsqtrue)*90/pi) + 1.0) ! "true" theta_pq
csk	  K = ((acos(cosPf)*90/pi) + 1.0)
	  if(K.gt.0 .and. K.le. 90) then
csk	  if(K.gt.0 .and. K.le. 90 .and. iselast) then
	    thpqdis(K) = thpqdis(K) + Weight
	    if(Cut1) thdis1cut(K) = thdis1cut(K) + Weight!*eldelphi
	    if(Cut2) thdis2cut(K) = thdis2cut(K) + Weight
	    if(Cut3) thdis3cut(K) = thdis3cut(K) + Weight!*SigP
	  endif
c	  print *, ' Fourth'
	  K = (Wstar/0.01 + 1.0)
	  if(K.gt.0 .and. K.le. 450) then
	    Wstardis(K,Ihel) = Wstardis(K,Ihel) + Weight
	    if(Cut1) Wsdis1cut(K,Ihel) = Wsdis1cut(K,Ihel) + Weight!*eldelphi
	    if(Cut2) Wsdis2cut(K,Ihel) = Wsdis2cut(K,Ihel) + Weight
	    if(Cut3) Wsdis3cut(K,Ihel) = Wsdis3cut(K,Ihel) + Weight!*SigP
	  endif
c	  print *, ' Fifth'
	  K = (Wlab/0.01 + 1.0)
	  if(K.gt.0 .and. K.le. 450) then
	    Wlabdis(K,Ihel) = Wlabdis(K,Ihel) + Weight
	    if(Cut1) Wldis1cut(K,Ihel) = Wldis1cut(K,Ihel) + Weight!*eldelphi
	    if(Cut2) Wldis2cut(K,Ihel) = Wldis2cut(K,Ihel) + Weight
	    if(Cut3) Wldis3cut(K,Ihel) = Wldis3cut(K,Ihel) + Weight!*SigP
	  endif
c	  print *, ' Sixth'

	  K = (Wlab/0.02 + 1.0)
	  if((K.gt.0) .and. (K.le. 450).and.Cut3) then !CHOOSE WHICH CUT YOU WANT HERE
	     L = 22+13.0*(alog(Q2/0.91883882)/alog(10.0)) ! EG1b Q2 bins
	    if (L.gt.0 .and. L.le.45) then
	      WQ2Sumdis(K,L) = WQ2Sumdis(K,L) + 1.0*Weight!*eldelphi
	      WQ2Diffdis(K,L) = WQ2Diffdis(K,L) + (3.0-2.0*Ihel)*Weight!*eldelphi
	    endif
	  endif
c	  print *, ' Made it through standard histos'
csk
CSK: CHANGED TO ALLOW ONLY DIS EVENTS and for finer bins at lo x

c	  if (xBjLab .ge. 0.1) then
c	    K = (xBjLab/0.05 + 1.0)
c	  else if (xBjlab .gt. 0.08) then
c	    K = 2
c	  else if (xBjlab .gt. 0.06) then
c	    K = 1
c	  else
c	    K = 0
c	  endif
ccc	  L = (alog(Q2/0.00681292)/alog(10.0))*3.0+1.0 ! (Q2/1.0 + 1.0)
c	  L = (alog(Q2)/alog(10.0))*3.0+7.0 ! (Q2/1.0 + 1.0)
c	  if((K.gt.0).and.(K.le.20).and.(L.gt.0).and.(L.le.15).and.
c     *		(Wlab .ge. 2.00).and.Cut1) then	! WLab
c	    xQ2dis(K,L) = xQ2dis(K,L) + Weight*eldelphi
c	    xQ2dis2(K,L) = xQ2dis2(K,L) + (3.0-2.0*Ihel)*Weight*eldelphi
c	  endif
csk	  
	  K = (xilab/0.01 + 1.0) ! Deuteron Hardcode !
c	  L = (alphamiss/0.1 + 1.0) ! measured alpha
	  L = (alphastrue/0.1 + 1.0) ! "true" alpha
	  if((K.gt.0 .and. K.le. 100).and.(L.gt.0 .and. L.le.20).and.Cut1)
     *		xialpdis(K,L) = xialpdis(K,L) + Weight!*eldelphi
	  K = (xinrs/0.01 + 1.0)
	  L = ((Mstar**2)/0.045 + 1.0)
	  if((K.gt.0 .and. K.le. 100).and.(L.gt.0 .and. L.le.20).and.Cut1)
     *		ximstdis(K,L) = ximstdis(K,L) + Weight!*eldelphi
	      
	else ! NOT I = 1
	  write (3,1999)
1999	  format(' ',/)

1001      format(f12.3,4e12.4)
	  write (3,1002)
1002	  format ('  Nu_L_[GeV] sigma_[cm2] sCut1_[cm2] sCut2_[cm2] sCut3_[cm2]')
          do K = 1 , 100
            write (3,1001) (K-0.5)*0.1, nudis(K),  
     *		nudis1cut(K), nudis2cut(K), nudis3cut(K)
          enddo
	  
2001      format(f12.3,8e12.4)
	  write (3,1999)
	  write (3,1003)
1003	  format ('  Q2_[GeV2]  sigma+[cm2] sigma-[cm2] sCut1+[cm2] sCut1-[cm2] ',
     &     'sCut2+[cm2] sCut2-[cm2] sCut3+[cm2] sCut3-[cm2]')
          do K = 1 , 150
            write (3,2001) (K-0.5)*0.1, Q2dis(K,1), Q2dis(K,2), Q2dis1cut(K,1), 
     *		Q2dis1cut(K,2), Q2dis2cut(K,1), Q2dis2cut(K,2), Q2dis3cut(K,1), Q2dis3cut(K,2)
          enddo
	  
csek2007 !!!!
csk	  goto 1847
csek2007 !!!!	  
	  
	  write (3,1999)
	  write (3,1004)
1004	  format ('  P_s_[GeV]  sigma_[cm2] sCut1_[cm2] sCut2_[cm2] sCut3_[cm2]')
csk1004      format ('  P_f_[GeV]  sigma_[cm2] sCut1_[cm2] sCut2_[cm2] sCut3_[cm2]')
          do K = 1 , 100
            write (3,1001) (K-0.5)*0.01, Psdis(K), 
     *		Psdis1cut(K), Psdis2cut(K), Psdis3cut(K)
csk            write (3,1001) (K-0.5)*0.1, Psdis(K), 
csk     *		Psdis1cut(K), Psdis2cut(K), Psdis3cut(K)
          enddo
	  
	  write (3,1999)
	  write (3,1005)
1005	  format ('  theta_qp_s sigma_[cm2] sCut1_[cm2] sCut2_[cm2] sCut3_[cm2]')
csk1005	  format ('  theta_qp_f sigma [cm2] sCut1_[cm2] sCut2_[cm2] sCut3_[cm2]')
          do K = 1 , 90
            write (3,1001) (K-0.5)*2.0, thpqdis(K),
     *		thdis1cut(K), thdis2cut(K), thdis3cut(K)
          enddo
	  
	  write (3,1999)
	  write (3,1006)
1006	  format ('  Wstar[GeV] sigm+_[cm2] +Cut1_[cm2] +Cut2_[cm2] +Cut3_[cm2]',
     *		' sigm-_[cm2] -Cut1_[cm2] -Cut2_[cm2] -Cut3_[cm2]')
          do K = 1 , 450
            write (3,1080) ((K-0.5)*0.01), (Wstardis(K,L),
     *		Wsdis1cut(K,L), Wsdis2cut(K,L), Wsdis3cut(K,L) , L = 1,2)
          enddo
	  
	  write (3,1999)
	  write (3,1007)
1007	  format ('  Wlab_[GeV] sigm+_[cm2] +Cut1_[cm2] +Cut2_[cm2] +Cut3_[cm2]',
     *		' sigm-_[cm2] -Cut1_[cm2] -Cut2_[cm2] -Cut3_[cm2]')
          do K = 1 , 450
            write (3,1080) ((K-0.5)*0.01), (Wlabdis(K,L),
     *		Wldis1cut(K,L), Wldis2cut(K,L), Wldis3cut(K,L) , L = 1,2)
          enddo
1080      format(f12.3,8e12.4)

	  write (3,1999)
	  write (3, 1019) ((L*0.5 - 0.25), L=1,30)
1019	  format(' Nu vs. Q2: ',/,
     *		'     Nu_/_Q2',30f12.4)
	  do K = 1 ,100
	    write(3,1030) ((K-0.5)*0.1),(nuQ2dis(K,L), L=1,30)
	  enddo
	  
c	  write(3,1999)
c	  write(3,901) ((float(Ip)/10.- 0.1), Ip = 1, 5)
c901	  format(' Elastic ps_vs_cos (Sum): ',/,'  IQ2          Q2    cos_/_ps', 5f12.4)
c	  do L = 16,34
c	  Q2 = 0.91883882*10.0**((float(L)-26.5)/13.0)
c	    do Icos = 1,3
c	      cosPmissq = 0.675*(Icos-2)
c	      write(3,909) L,Q2, cosPmissq,(mikeDis(L,Ip,Icos,1), Ip = 1,5)
c	    enddo
c	  enddo
c909	  format(' ',I4, 2f12.4,5e12.4)
c	  write(3,902) ((float(Ip)/10.- 0.1), Ip = 1, 5)
c902	  format(' Elastic ps_vs_cos (Diff): ',/,'  IQ2          Q2    cos_/_ps', 5f12.4)
c	  do L = 16,34
c	    Q2 = 0.91883882*10.0**((float(L)-26.5)/13.0)
c	    do Icos = 1,3
c	      cosPmissq = 0.675*(Icos-2)
c	      write(3,909) L, Q2, cosPmissq,(mikeDis(L,Ip,Icos,2), Ip = 1,5)
c 	    enddo
c	  enddo
	  
	  write (3,1999)
	  write (3, 1021) ( (10.0**((float(L)-22.0)/13.0)), L=21,44)
1021	format(' Wlab_vs_Q2 (Sum): ',/,' W_/_Q2 ',24f12.4)
	    do K = 1 ,225
	      write(3,1030) ((K-0.5)*0.02),(WQ2Sumdis(K,L), L=21,44)
	    enddo

	  write (3,1999)
	  write (3, 1022) ( (10.0**((float(L)-22.0)/13.0)), L=21,44)
1022	format(' Wlab_vs_Q2 (Diff): ',/,' W_/_Q2 ',24f12.4)
	    do K = 1 ,225
	      write(3,1030) ((K-0.5)*0.02),(WQ2Diffdis(K,L), L=21,44)
	    enddo

CSK: CHANGED TO ALLOW ONLY DIS EVENTS


1847      write (3,1999)
ccc          write(3, 1041) ( (0.01*10.0**((float(L)-1.0)/3.0)), L = 1,15)
          write(3, 1041) ( (10.0**((float(L)-6.5)/3.0)), L = 1,15)
1041      format(' xBjL_vs_Q2 sum: ',/,'x_/_Q2',15f12.4)
          write(3,1040) (0.07),(xQ2dis(1,L), L=1,15)
          write(3,1040) (0.09),(xQ2dis(2,L), L=1,15)
          do K = 3 ,20
            write(3,1040) ((K-0.5)*0.05),(xQ2dis(K,L), L=1,15)
          enddo

	  write (3,1999)
ccc          write(3, 1042) ( (0.01*10.0**((float(L)-1.0)/3.0)), L = 1,15)
          write(3, 1042) ( (10.0**((float(L)-6.5)/3.0)), L = 1,15)
1042      format(' xBjL_vs_Q2 diff: ',/,'x_/_Q2',15f12.4)
          write(3,1040) (0.07),(xQ2dis2(1,L), L=1,15)
          write(3,1040) (0.09),(xQ2dis2(2,L), L=1,15)
          do K = 3 ,20
            write(3,1040) ((K-0.5)*0.05),(xQ2dis2(K,L), L=1,15)
          enddo

	  write (3,1999) ! Deuteron Hardcode !
	  write (3, 1011) (((L-0.5)*0.1), L=1,20)
1011	  format(' xi_Nachtmann(lab) vs. alpha_S: ',/,
     *		'xi_/_alpha_S',20f12.4)
	  do K = 1 ,100
	    write(3,1010) ((K-0.5)*0.01),(xialpdis(K,L), L=1,20)
	  enddo

	  write (3,1999) ! Deuteron Hardcode !
	  write (3, 1012) ((sqrt((L-0.5)*0.045)), L=1,20)
1012	  format(' xi_Nachtmann(true) vs. MStar: ',/,
     *		'xi_n_/_MStar',20f12.4)
	  do K = 1 ,100
	    write(3,1010) ((K-0.5)*0.01),(ximstdis(K,L), L=1,20)
	  enddo
	  
1010	format(f12.4,20e12.4)
1040	format(f12.4,15e12.4)
1030	format(f12.4,44e12.4)	
	endif     ! I = 1                                                   
	return
	end
