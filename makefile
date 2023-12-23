#FLAGS =   -w 
#FLAGS = -c  -w 
#FLAGS =   -C -w -save
#FLAGS = -O3 -qextchk -qarch=auto 
#FLAGS = -O3  -qarch=auto 
#FLAGS = -g -d1 -C -w  -qattr 

 FC = xlf_r
#FC = gfortran


 STRUCT = groucho_gapbld.o  gettime.o  dexptablesmall_setup.o dexptablebig_setup.o synaptic_map_construct.o synaptic_compmap_construct.o durand.o 

 INTEGRATE =  integrate_LOT.o fnmda.o integrate_supVIP.o integrate_deepLTS.o integrate_deepbask.o integrate_deepng.o  integrate_placeholder6.o integrate_placeholder3.o integrate_placeholder2.o integrate_placeholder1.o integrate_supng.o   integrate_placeholder5.o    integrate_placeholder4.o 


piriformENDO : piriformENDO.f makefile gettime.o durand.o $(STRUCT) $(INTEGRATE)
	mpif90 -g -O3 -qfixed piriformENDO.f $(STRUCT) $(INTEGRATE) -o piriformENDO

gettime.o : gettime.c
	gcc -c -O2 gettime.c

durand.o : durand.f
	xlf_r -c -O3 durand.f

%.o : %.f
	$(FC)  -c -g -O3 $<


clean :
	rm -f piriformENDO    
	
