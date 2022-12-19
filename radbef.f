	subroutine radbefB12(I, Q2, Ein, Eout, depol)
	
C	Subroutine to apply external radiation loss to incoming electron
C	Created 10-Jul-2003 Sebastian Kuhn
C   Updated Dec-2022 to cludge ALL radiative effects - SK

	IMPLICIT NONE
	REAL*4 ran1, Ein, Eout, K, RanNum, temp, Q2, teqrad
csk		new variables for depolarization effect
	real*4 depol, PSI1 /18.2309/, PSI2 /17.5642/
csk
	REAL*4 Kmin /0.001/	! Ignore radiative losses below 0.001 GeV
c       SK: The following parameters can be gotten from the Excel spreadsheet PreRad_eloss.xls
c	REAL*4 PKmin /0.00551/	! Probability of loss>Kmin (for 4.2 GeV BoNuS)
c	REAL*4 const1 /0.00775/, const2 /0.092/	! radiation loss parameters, product = approx. bt
	INTEGER I
	REAL*4 PKmin	! Probability of loss>Kmin (for 2 GeV EG4)
	REAL*4 const1, const2 /0.05245/	! radiation loss parameters, product = approx. bt
    	REAL*4 aoffset /0.1266/, aslope /16.2249/, pkminoffset /0.0564/, pkminslope /7.6466/
    
c    	write(6,1) Ein, Q2
1   	format(' BEFORE:  Ein: ',F9.5,' Q2,',F9.5)
    	Eout = Ein
    	if(Ein.lt.1.0) return

    	teqrad = (3.0/4.0)*(1.0/137.036/3.1415927)*(alog(Q2/0.000511/0.000511)-1.0)
    	const1 = aoffset + aslope*teqrad
    	PKmin = pkminoffset + pkminslope*teqrad
c    	write(6,2) teqrad, const1, PKmin
2   	format(' Eq.Rad. t: ',f9.5,' a-parameter: ',f9.5,' Prob. for kmin: ',f9.5)

123	RanNum = ran1(I)
c	write(6,5) RanNum
5	format(' Random Number: ',f9.5)
	if (RanNum .ge. PKmin) return
	
	temp = (PKmin - RanNum)/const1
c	write(6,6) temp
6	format(' first temp:',f9.5)
	if ((temp .lt. 0.0) .or. (temp .gt. 0.9999)) goto 123 ! return
	
	temp = asin(temp)/const2
c	write(6,7) temp
7	format(' second temp: ',f9.5)
	K = Kmin*exp(temp)
	Eout = Ein - K
c    	write(6,3) K, Eout
3   	format(' Photon energy:',f9.5,' Eout:',f9.5)

	if (Eout .le. 0.0) goto 123 !Eout = 0.001
csk	put in depolarization
	depol = K*K*(PSI2*2./3.) /
     *		((Ein*Ein + Eout*Eout)*PSI1 - Ein*Eout*PSI2*2./3.)
c    	write(6,4) depol
4   	format(' Depolarization: ',f9.5)

	return
	end
	
