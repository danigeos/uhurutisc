#rm /home/ivone/MAIN/src/Uhuru/src/*.o /home/ivone/MAIN/src/tisc/src/call_tisc.o /home/ivone/MAIN/src/tisc/src/surf_proc.o /home/ivone/MAIN/src/tisc/lib/libreria.o

FC=gfortran   ## Intel	## /tisc/src/call_tisc_ifort.c :	
			## 	int init_tisc_(
			##	int call_surf_proc_(
			##	int call_flexure_(
### -------------------------------------------------------------------------------------------------

#FC=f90   ## Absoft	## /tisc/src/call_tisc_f90.c :	
			## 	int INIT_TISC (
			##	int CALL_SURF_PROC (
			##	int CALL_FLEXURE (
### -------------------------------------------------------------------------------------------------

DIR_main=/Users/danielgc/software/uhuru/Uhuru
#DIR_main=`pwd`

BIN_UHU = UHURU
BINDIR = bin
SRCDIR = src
LIBMATH = lib/mathematics.f
LIBOUTIN = lib/outin.f
LIBYIELD = lib/Ch_tree_lib.f
LIBT1D = lib/Termica_1D_lib.f
LIBTH_RW = lib/RW_PARAMETRES.f
LIBTH_noe = lib/RW_PARAMETRES_noelev.f

DIR_tisc = $(DIR_main)/tisc
LIBsurf = $(DIR_tisc)/src/surf_proc.o
LIBtisc = $(DIR_tisc)/lib/libreria.o
LIBcalltisc = $(DIR_tisc)/src/call_tisc.o

all: thin_sheet.o UHURU.o $(BINDIR)/$(BIN_UHU) $(BINDIR)/graficsth

####################################################################################
$(BINDIR)/graficsth: $(SRCDIR)/graficsth.f $(LIBTH_RW)
	$(FC) -o $(BINDIR)/graficsth $(SRCDIR)/graficsth.f $(LIBTH_RW) $(LIBOUTIN)

call_tisc.o: 
	(cd $(DIR_tisc)/src; make; make call_tisc.o; make surf_proc.o; cd $(DIR_tisc)/lib; make libreria.o)

thin_sheet.o: $(SRCDIR)/thin_sheet.f $(LIBMATH)
	$(FC) -c $(SRCDIR)/thin_sheet.f $(LIBMATH)

UHURU.o: $(SRCDIR)/UHURU.f $(LIBTH_RW) $(LIBOUTIN) $(LIBYIELD) $(LIBT1D)
	$(FC) -c $(SRCDIR)/UHURU.f $(LIBTH_RW) $(LIBOUTIN) $(LIBYIELD) $(LIBT1D)	

$(BINDIR)/$(BIN_UHU): call_tisc.o thin_sheet.o UHURU.o
	$(FC) -lc -o $(BINDIR)/$(BIN_UHU) $(LIBsurf) $(LIBcalltisc) $(LIBtisc) *.o -lm



