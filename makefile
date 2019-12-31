#FFLAGS = -g -fast -ftrap=common
#FFLAGS = -g -ftrap=common
#FFLAGS =  -O vector0 -O scalar3 -O inline3
#FFLAGS = -Ofast 
FFLAGS = -O3
FC = gfortran

selfcons : sor.o selfcons.o setupbcphi.o init.o jangd.o evalrho.o getstripmass.o getcloudmass.o getphic.o getthermener.o getbindener.o getrotener.o getvolume.o getvphi.o dump.o
	gfortran $(FFLAGS) -o selfcons sor.o selfcons.o setupbcphi.o init.o jangd.o evalrho.o getstripmass.o getcloudmass.o getphic.o getthermener.o getbindener.o getrotener.o getvolume.o getvphi.o dump.o

#.F.o:
#	gcc -E -P $< > $*.f
#	f77 ${FFLAGS} -c $*.f
 
#----------------------------------------------------------------------------
clean :
	rm -f *.o *.f
	rm -f selfcons
	rm -f core
#----------------------------------------------------------------------------
newrun :
	rm -f isorot* fort.*
