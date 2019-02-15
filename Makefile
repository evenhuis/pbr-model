#_a makefile



#FC= gfortran-mp-4.5 -O2   -Warray-bounds -O0 -g -fbounds-check
FC= gfortran -std=f2008  -O2  # -Warray-bounds -p -O0 -g -fbounds-check

INC=/Users/evenhuis/Dropbox/fortran_programs/lib/
BLAS= -L/usr/lib -lblas -llapack

LIB_DIR = /Users/evenhuis/Dropbox/fortran_programs/lib/
LIB_O   = -I$(LIB_DIR) -L$(LIB_DIR) -lmylib 

%: %.f90 util_mod.o
	$(FC) $(BLAS) -o $* $*.f90   $(LIB_O) util_mod.o

# - - - - - - - growth_atten 1 - - - -  - --  - - --  - - - -
growth_att_mod.o : growth_att_mod.f90                                           util_mod.o ext_driver_mod.o
	$(FC) $(BLAS) -c                                         growth_att_mod.f90                           $(LIB_O)
	$(FC) $(BLAS) -c                       interface_mod.f90 growth_att_mod.o util_mod.o ext_driver_mod.o $(LIB_O)            
	$(FC) $(BLAS) -c  py_interface_mod.f90 interface_mod.o   growth_att_mod.o util_mod.o ext_driver_mod.o $(LIB_O)

growth_att:                     growth_att.f90                      growth_att_mod.o util_mod.o ext_driver_mod.o
	$(FC) $(BLAS) -o growth_att  growth_att.f90 interface_mod.o   growth_att_mod.o util_mod.o ext_driver_mod.o $(LIB_O)

growth_atten.so: growth_att_mod.o
	touch growth_att.so
	rm    growth_att.so
	f2py -c -m growth_atten py_interface_mod.f90 growth_att_mod.o interface_mod.o ext_driver_mod.o util_mod.o $(LIB_O) $(BLAS)


# - - - - - - - Model 1 - - - -  - --  - - --  - - - -
pr_model1_mod.o : pr_model1_mod.f90  														 util_mod.o ext_driver_mod.o
	$(FC) $(BLAS) -c            			    						pr_model1_mod.f90 								  $(LIB_O)
	$(FC) $(BLAS) -c 							   interface_mod.f90 pr_model1_mod.o util_mod.o ext_driver_mod.o $(LIB_O)            
	$(FC) $(BLAS) -c  py_interface_mod.f90 interface_mod.o   pr_model1_mod.o util_mod.o ext_driver_mod.o $(LIB_O)

pr_model1:                     pr_model1.f90                   pr_model1_mod.o util_mod.o ext_driver_mod.o
	$(FC) $(BLAS) -o pr_model1  pr_model1.f90 interface_mod.o   pr_model1_mod.o util_mod.o ext_driver_mod.o $(LIB_O)

pr_model1.so: pr_model1_mod.o
	touch pr_model1.so
	rm    pr_model1.so
	f2py -c -m pr_model1 py_interface_mod.f90 pr_model1_mod.o interface_mod.o ext_driver_mod.o util_mod.o $(LIB_O) $(BLAS)

# - - - - - - - Model 2 - - - -  - --  - - --  - - - -
pr_model2_mod.o : pr_model2_mod.f90  														 util_mod.o ext_driver_mod.o
	$(FC) $(BLAS) -c            			    						pr_model2_mod.f90 								  $(LIB_O)
	$(FC) $(BLAS) -c 							   interface_mod.f90 pr_model2_mod.o util_mod.o ext_driver_mod.o $(LIB_O)            
	$(FC) $(BLAS) -c  py_interface_mod.f90 interface_mod.o   pr_model2_mod.o util_mod.o ext_driver_mod.o $(LIB_O)

pr_model2:                     pr_model2.f90                   pr_model2_mod.o util_mod.o ext_driver_mod.o
	$(FC) $(BLAS) -o pr_model2  pr_model2.f90 interface_mod.o   pr_model2_mod.o util_mod.o ext_driver_mod.o $(LIB_O)

pr_model2.so: pr_model2_mod.o
	touch pr_model2.so
	rm    pr_model2.so
	f2py -c -m pr_model2 py_interface_mod.f90 pr_model2_mod.o interface_mod.o ext_driver_mod.o util_mod.o $(LIB_O) $(BLAS)

# - - - - - - - Model 6 - - - -  - --  - - --  - - - -
pr_model6_mod.o : pr_model6_mod.f90  							      ext_driver_mod.o
	$(FC) $(BLAS) -c               			   pr_model6_mod.f90 ext_driver_mod.o $(LIB_O)
	$(FC) $(BLAS) -c  py_pr_interface_mod.f90 pr_model6_mod.o   ext_driver_mod.o $(LIB_O)

pr_model6.so: pr_model6_mod.o
	touch pr_model6.so
	rm    pr_model6.so
	f2py -c -m pr_model6 py_pr_interface_mod.f90 pr_model6_mod.o ext_driver_mod.o  $(LIB_O) $(BLAS)

pr_model6:                     pr_model6.f90                   pr_model6_mod.o util_mod.o ext_driver_mod.o
	$(FC) $(BLAS) -o pr_model6  pr_model6.f90                   pr_model6_mod.o util_mod.o ext_driver_mod.o $(LIB_O)



# - - - - - - - Model 8 - - - -  - --  - - --  - - - -
pr_model8_mod.o : pr_model8_mod.f90  							      ext_driver_mod.o
	$(FC) $(BLAS) -c               			   pr_model8_mod.f90 ext_driver_mod.o $(LIB_O)
	$(FC) $(BLAS) -c  py_pr_interface_mod.f90 pr_model8_mod.o   ext_driver_mod.o $(LIB_O)

pr_model8.so: pr_model8_mod.o 
	touch pr_model8.so
	rm    pr_model8.so
	f2py -c -m pr_model8 py_pr_interface_mod.f90 pr_model8_mod.o ext_driver_mod.o  $(LIB_O) $(BLAS)

# - - - - - - - Model 9 - - - -  - --  - - --  - - - -
pr_model9_mod.o : pr_model9_mod.f90                            ext_driver_mod.o
	$(FC) $(BLAS) -c                           pr_model9_mod.f90 ext_driver_mod.o $(LIB_O)
	$(FC) $(BLAS) -c  py_pr9_interface_mod.f90 pr_model9_mod.o   ext_driver_mod.o $(LIB_O)

pr_model9.so: pr_model9_mod.o
	touch pr_model9.so
	rm    pr_model9.so
	f2py -c -m pr_model9 py_pr9_interface_mod.f90 pr_model9_mod.o ext_driver_mod.o  $(LIB_O) $(BLAS)


# - - - - - - - Model 10 - - - -  - --  - - --  - - - -
pr_model10_mod.o : pr_model10_mod.f90                            ext_driver_mod.o
	$(FC) $(BLAS) -c                           pr_model10_mod.f90 ext_driver_mod.o $(LIB_O)
	$(FC) $(BLAS) -c  py_pr_model10_mod.f90 pr_model10_mod.o   ext_driver_mod.o $(LIB_O)

pr_model10.so: pr_model10_mod.o
	touch pr_model10.so
	rm    pr_model10.so
	f2py -c -m pr_model10 py_pr_model10_mod.f90 pr_model10_mod.o ext_driver_mod.o  $(LIB_O) $(BLAS)

# - - - - - - - Biomass - - - -  - --  - - --  - - - -
biomass_model_mod.o : biomass_model_mod.f90  							      ext_driver_mod.o
	$(FC) $(BLAS) -c               			   biomass_model_mod.f90 ext_driver_mod.o $(LIB_O)

biomass_model.so: biomass_model_mod.o 
	touch biomass_model.so
	rm    biomass_model.so
	$(FC) $(BLAS) -c     py_biomass_mod.f90 biomass_model_mod.o ext_driver_mod.o $(LIB_O)
	f2py -c -m biomass_model py_biomass_mod.f90 biomass_model_mod.o ext_driver_mod.o $(LIB_O) $(BLAS)

biomass_model:                     biomass_model.f90     biomass_model_mod.o  ext_driver_mod.o
	$(FC) $(BLAS) -o biomass_model  biomass_model.f90     biomass_model_mod.o  ext_driver_mod.o $(LIB_O)



test_intp: test_intp.f90 util_mod.o
	$(FC) $(BLAS) -o test_intp test_intp.f90 util_mod.o  $(LIB_O)

test_light: test_light.f90 ext_driver_mod.o
	$(FC) $(BLAS) -o test_intp test_intp.f90 ext_driver_mod.o  $(LIB_O)

util_mod.o: util_mod.f90
	$(FC) $(BLAS) -c util_mod.f90 $(LIB_O)

util_new_mod.o: util_new_mod.f90
	$(FC) $(BLAS) -c util_new_mod.f90 $(LIB_O)

ext_driver_mod.o : ext_driver_mod.f90
	 $(FC) $(BLAS) -c ext_driver_mod.f90 $(LIB_O)





