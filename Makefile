FFLAGS = -c -ffixed-line-length-none -fno-automatic \
         -O2 -fno-f2c -funroll-loops \
         -Wimplicit -Wsurprising 

# not used:   -fno-silent -pedantic -Wall

# NOTE: use of platform optimizing flags (such as -m486) is STRONGLY
# recommended -- check gcc (for g77) info or man pages

# description of flags:
# fno-silent, pedantic, Wimplicit, Wall and Wsurprising  generate warnings
#  for bad coding -- just don't do it!
# fixed-line-length-none  precludes source code lines being truncated
# O2  causes compiler to take time to optimize for speedy running
# malign-double  aligns variables with word boundaries on PCs, is faster
# m486  optimize code for i486 platform (as opposed to i386) -- no Pentium yet
# funroll-loops  linearizes short loops, avoiding the branching & counting
#  delays at runtime
# c  compile only, do not link
# fno-f2c  compiled objects are f77 INCOMPATIBLE -- only an issue if linking
#  to objects or libraries compiled with f77
# fno-automatic  makes all variables static (SAVE) instead of undefining
#  them when exiting the subroutine etc.
#  NOTE: this is needed until the respective variables are SAVEd individually
#        because the XLF compiler in AIX does it and that's where the code
#        was developed  the flags used there were:  
#        -fast -O2  -c -qdpc -qautodbl
#        the meaning of qdpc is to turn all reals into real*8
#        qautodbl allows singles to be computed as doubles (speed)


FC = gfortran

LIB = pol.a
INCLUDE = radcon.inc radvarfix.inc radvar.inc prinplot.inc instruct.inc

OALL =	${LIB}(h2model_thia.o) ${LIB}(models.o) ${LIB}(i_d2_model.o) ${LIB}(sek_me.o) \
	${LIB}(r1998.o) ${LIB}(sek_nce.o) ${LIB}(sek_fit.o) ${LIB}(frw_fit.o) ${LIB}(sek_2021.o) \
	${LIB}(RRicco.o) ${LIB}(crosssection.o) ${LIB}(ran1.o) ${LIB}(f1f221redux.o) \
	${LIB}(cross_radel.o) ${LIB}(cross_radin.o) ${LIB}(elasradcor.o) \
	${LIB}(movingcross.o) ${LIB}(writeout.o) ${LIB}(writeouteg4.o) \
	${LIB}(mocross_radel.o) ${LIB}(mocross_radin.o) ${LIB}(cross_radeut.o) \
	${LIB}(protontrack.o) ${LIB}(radbef.o) ${LIB}(tailread.o) ${LIB}(ratioread.o) \
	${LIB}(elasradtail.o) ${LIB}(eqrad.o) ${LIB}(clas12acc.o) \
	${LIB}(F1F209.o) ${LIB}(F1F2IN06.o) ${LIB}(F1F2QE06.o) \
	${LIB}(f1f209_dr.o) ${LIB}(rdcorr.o) ${LIB}(F1F207.o) ${LIB}(lithium.o) \
	${LIB}(sf.o) ${LIB}(resmod.o) ${LIB}(H2HallC.o) ${LIB}(christy.o) \
	${LIB}(f2allm.o) ${LIB}(radaft.o) ${LIB}(newSFs.o) ${LIB}(DSFs.o) \
	${LIB}(deuteron_new.o) ${LIB}(lcsupp_new.o) ${LIB}(Yoni.o) ${LIB}(g1Dfit.o) \
	${LIB}(movingcross_MEIC.o) ${LIB}(writeout_MEIC.o) ${LIB}(writeout_clas12.o)
	
	
#	${LIB}(F1F207.o) ${LIB}(f2allm.o) ${LIB}(radaft.o) ${LIB}(newSFs.o) \
#	Choose one of the two versions only!
$(OALL): $(INCLUDE)
maindasym.o : $(INCLUDE)
dgampn.o : $(INCLUDE)
maindasym_meic.o : $(INCLUDE)
maindasym_clas12.o : $(INCLUDE)
mainLasym_clas12.o : $(INCLUDE)
sigma.o : $(INCLUDE)
sigma_el.o : $(INCLUDE)
strucfunc.o : $(INCLUDE)
wf.o : $(INCLUDE)
test.o : $(INCLUDE)
radcheck.o : $(INCLUDE)
main_elasrad.o : $(INCLUDE)
main_elasintrad.o : $(INCLUDE)
main_radin.o : $(INCLUDE)
main_radtest.o : $(INCLUDE)
main_radeut.o : $(INCLUDE)
maindasym_cc.o : $(INCLUDE)
maindasym_f2ps.o : $(INCLUDE)
maindasym_f2s.o : $(INCLUDE)
SFextract.o : $(INCLUDE)
SFextract_cc.o : $(INCLUDE)
strucfunc.o : $(INCLUDE)
strucfunc_wm.o : $(INCLUDE)
lowxx.o : $(INCLUDE)
A1test.o : $(INCLUDE)
Moment.o : $(INCLUDE)
elasasym.o : $(INCLUDE)
testdwf.o : $(INCLUDE)
main_simplerad.o : $(INCLUDE)
pol.a: ${OALL}

deut: pol.a maindasym.o
	$(FC) -g -o deut.exe maindasym.o pol.a

dgampn: pol.a dgampn.o
	$(FC) -g -o dgampn.exe dgampn.o pol.a

deutmeic: pol.a maindasym_meic.o
	$(FC) -g -o deutmeic.exe maindasym_meic.o pol.a

deut12: pol.a maindasym_clas12.o
	$(FC) -g -o clas12.exe maindasym_clas12.o pol.a

lith12: pol.a mainLasym_clas12.o
	$(FC) -g -o polEMClm.exe mainLasym_clas12.o pol.a

deut_f2ps: pol.a maindasym_f2ps.o
	$(FC) -g -o deut_f2ps.exe maindasym_f2ps.o pol.a

deut_f2: pol.a maindasym_f2.o
	$(FC) -g -o deut_f2.exe maindasym_f2.o pol.a

deut_cc: pol.a maindasym_cc.o
	$(FC) -g -o deut_cc.exe maindasym_cc.o pol.a

sigma: pol.a sigma.o
	$(FC) -g -o sigma.exe sigma.o pol.a
	
sigma_el: pol.a sigma_el.o
	$(FC) -g -o sigma_el.exe sigma_el.o pol.a
	
wf:	pol.a wf.o
	$(FC) -g -o wf.exe wf.o pol.a

test: 	pol.a test.o
	$(FC) -g -o test.exe test.o pol.a

radin: 	pol.a main_radin.o
	$(FC) -g -o radin.exe main_radin.o pol.a

radtest: 	pol.a main_radtest.o
	$(FC) -g -o radtest.exe main_radtest.o pol.a

radeut: 	pol.a main_radeut.o
	$(FC) -g -o radeut.exe main_radeut.o pol.a

radcheck: 	pol.a radcheck.o
	$(FC) -g -o radcheck.exe radcheck.o pol.a
	
elasrad:	pol.a main_elasrad.o
	$(FC) -g -o elasrad.exe main_elasrad.o pol.a
	
elasintrad:	pol.a main_elasintrad.o
	$(FC) -g -o elasintrad.exe main_elasintrad.o pol.a
	
SF: 		pol.a SFextract.o
	$(FC) -g -o SF.exe SFextract.o pol.a
	
SFcc: 		pol.a SFextract_cc.o
	$(FC) -g -o SFcc.exe SFextract_cc.o pol.a
	
strucfunc: 	pol.a strucfunc.o
	$(FC) -g -o strucfunc.exe strucfunc.o pol.a
	
strucfuncwm: 	pol.a strucfunc_wm.o
	$(FC) -g -o strucfunc_wm.exe strucfunc_wm.o pol.a
	
lowxx: 	pol.a lowxx.o
	$(FC) -g -o lowxx.exe lowxx.o pol.a
	
A1test: 	pol.a A1test.o
	$(FC) -g -o A1test.exe A1test.o pol.a
	
Moment: 	pol.a Moment.o
	$(FC) -g -o Moment.exe Moment.o pol.a
	
testFF: 	pol.a elasasym.o
	$(FC) -g -o testFF.exe elasasym.o pol.a

testdwf: 	pol.a testdwf.o
	$(FC) -g -o testdwf.exe testdwf.o pol.a
	
simplerad:	pol.a main_simplerad.o
	$(FC) -g -o main_simplerad.exe main_simplerad.o pol.a

