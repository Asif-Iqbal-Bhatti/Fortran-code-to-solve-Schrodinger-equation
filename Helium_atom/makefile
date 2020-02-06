# FMKMF_F90 not set: using f90
# FMKMF_SFTAG not set: using f90
# FMKMF_SPATH not set: using . 
# FMKMF_LINKOPTS not set: using no link options 
# Main program is HF.f90 
# process_fsource called with arg HF.f90 
# HF.f90 Uses Module matrixmultiplication
# Full list of modules in HF.f90: matrixmultiplication 
# Uses matrixmultiplication which is in ./matrixmulti.f90
# process_fsource called with arg ./matrixmulti.f90 
# Full list of modules in ./matrixmulti.f90:  

# ------------------Macro-Defs---------------------
F90=ifort -mkl=parallel 

# -------------------End-macro-Defs---------------------------

# Here is the link step 
HF:matrixmulti.o HF.o 
	 $(F90) -o HF matrixmulti.o HF.o   

# Here are the compile steps
 
matrixmulti.o:./matrixmulti.f90  
	 $(F90) -c ./matrixmulti.f90 

HF.o:HF.f90 matrixmulti.o 
	 $(F90) -c HF.f90 
# This entry allows you to type " make clean " to get rid of
# all object and module files 
clean:
	rm -f -r f_{files,modd}* *.o *.mod *.M *.d V*.inc *.vo \
	V*.f *.dbg album F.err
  
