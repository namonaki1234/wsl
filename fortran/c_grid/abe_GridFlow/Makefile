########################################################################################
######## マクロ定義								########
########################################################################################
FC = ifort
#FC = gfortran
OPT_O = -o
OPT_C = -c
FLAGS = -O3 -qopenmp  -mcmodel=large -traceback
# -traceback -check bounds
#FLAGS = -fbounds-check -ffpe-trap=invalid,zero,overflow -fbacktrace
#FLAGS = -O2
#FLAGS = -fast
#FLAGS = -CB -fpe0 -traceback -check -g
#FLAGS = -mcmodel=large -shared-intel
#FLAGS = -fast -openmp
FLAGS2 = -heap-arrays -qopenmp -O3
SRCDIR_GRID = ./src/GridFlow
SRCDIR_MPSINI = ./src/MPSini
SRCDIR_MPS = ./src/MPS
SRCDIR_ICM = ./src/ICM
MODDIR_GRID = ./mod/GridFlow
MODDIR_MPSINI = ./mod/MPSini
MODDIR_MPS = ./mod/MPS
OBJDIR_GRID = ./temp/GridFlow
OBJDIR_MPSINI = ./temp/MPSini
OBJDIR_MPS = ./temp/MPS
OBJDIR_ICM = ./temp/ICM
EXE = exe
OBJ = o
DEL = rm
#DEL = del
ECHO   = echo
ALLEXE = $(OBJDIR_GRID)/Grid.$(EXE) \
         $(OBJDIR_GRID)/Overset.$(EXE) \
         $(OBJDIR_GRID)/IniFlow.$(EXE) \
         $(OBJDIR_GRID)/Flow.$(EXE) \
         $(OBJDIR_GRID)/ViewFlow.$(EXE) \
         $(OBJDIR_MPSINI)/MPS_ini.$(EXE) \
         $(OBJDIR_MPS)/MPS.$(EXE) \
	 $(OBJDIR_ICM)/ICM.$(EXE)
GRIDEXE = $(OBJDIR_GRID)/Grid.$(EXE) \
          $(OBJDIR_GRID)/Overset.$(EXE) \
          $(OBJDIR_GRID)/IniFlow.$(EXE) \
          $(OBJDIR_GRID)/Flow.$(EXE) \
          $(OBJDIR_GRID)/ViewFlow.$(EXE)
MPSINIEXE = $(OBJDIR_MPSINI)/MPS_ini.$(EXE)
MPSEXE = $(OBJDIR_MPS)/MPS.$(EXE)
ICMEXE = $(OBJDIR_ICM)/ICM.$(EXE)

########################################################################################
########## サフィックスルール							########
########################################################################################
#.SUFFIXES: $(OBJ) .f90
#.f90.o:
#	$(FC) $(FLAGS) $(OPT_C) $<
# サフィックスルールは使えなさそう?

########################################################################################
########## 依存関係とコマンド							########
########################################################################################
# デフォルトターゲット
.PHONY: all
all: clean $(ALLEXE)

# 格子法
.PHONY: grid
grid: $(GRIDEXE)

# 粒子法initial
.PHONY: mpsini
mpsini:	$(MPSINIEXE)

# 粒子法
.PHONY: mps
mps: $(MPSEXE)

# 整理用
.PHONY: clean
clean :
	$(DEL) -f $(OBJDIR_GRID)/*.$(OBJ) $(OBJDIR_MPSINI)/*.$(OBJ) $(OBJDIR_MPS)/*.$(OBJ)
	$(DEL) -f $(OBJDIR_GRID)/*.$(EXE) $(OBJDIR_MPSINI)/*.$(EXE) $(OBJDIR_MPS)/*.$(EXE)
	$(DEL) -f *.mod

# 実行ファイル
$(OBJDIR_GRID)/Grid.$(EXE)	:	$(OBJDIR_GRID)/Grid_Generation_NACA.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_FileIO.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_NACA.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_Equation.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_Grid.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^
$(OBJDIR_GRID)/Overset.$(EXE)	:	$(OBJDIR_GRID)/Overset_Interpolation_NACA.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_FileIO.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_NACA.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_Grid.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^
$(OBJDIR_GRID)/IniFlow.$(EXE)	:	$(OBJDIR_GRID)/Initial_Flow_NACA.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_FileIO.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_NACA.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_Grid.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_Flow.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_ViewFlow.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^
$(OBJDIR_GRID)/Flow.$(EXE)	:	$(OBJDIR_GRID)/Flow_Field_NACA.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_FileIO.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_NACA.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_Flow.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_Icing.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^
$(OBJDIR_GRID)/ViewFlow.$(EXE)	:	$(OBJDIR_GRID)/View_Flow_NACA.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_FileIO.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_NACA.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_Flow.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_ViewFlow.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^

$(OBJDIR_MPSINI)/MPS_ini.$(EXE)	:	$(OBJDIR_MPSINI)/airfoil_2D.$(OBJ) \
					$(OBJDIR_MPSINI)/mod_is_2D.$(OBJ) \
					$(OBJDIR_MPSINI)/mod_is.$(OBJ) \
					$(OBJDIR_MPSINI)/mod_is_airfoil_2D.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^

$(OBJDIR_MPS)/MPS.$(EXE)	:	$(OBJDIR_MPS)/mps_icing_2D.$(OBJ) \
					$(OBJDIR_MPS)/mod_mps_2D.$(OBJ) \
					$(OBJDIR_MPS)/mod_mps.$(OBJ) \
					$(OBJDIR_MPS)/mod_mps_icing_2D.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_FileIO.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_NACA.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_Flow.$(OBJ) \
					$(OBJDIR_GRID)/Mod_Package_Grid.$(OBJ) \
					$(OBJDIR_MPS)/mod_package_droplet.$(OBJ)

	$(FC) $(FLAGS) $(OPT_O) $@ $^

$(OBJDIR_ICM)/ICM.$(EXE)	:	$(OBJDIR_ICM)/DeleteInnerParticle.$(OBJ)
	$(FC) $(FLAGS2) $(OPT_O) $@ $^

# オブジェクト
$(OBJDIR_GRID)/Grid_Generation_NACA.$(OBJ)	:	$(SRCDIR_GRID)/Grid_Generation_NACA.f90 \
							$(OBJDIR_GRID)/Mod_Package_FileIO.$(OBJ) \
							$(OBJDIR_GRID)/Mod_Package_NACA.$(OBJ) \
							$(OBJDIR_GRID)/Mod_Package_Equation.$(OBJ) \
							$(OBJDIR_GRID)/Mod_Package_Grid.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<
$(OBJDIR_GRID)/Overset_Interpolation_NACA.$(OBJ):	$(SRCDIR_GRID)/Overset_Interpolation_NACA.f90 \
							$(OBJDIR_GRID)/Mod_Package_FileIO.$(OBJ) \
							$(OBJDIR_GRID)/Mod_Package_NACA.$(OBJ) \
							$(OBJDIR_GRID)/Mod_Package_Grid.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<
$(OBJDIR_GRID)/Initial_Flow_NACA.$(OBJ)		:	$(SRCDIR_GRID)/Initial_Flow_NACA.f90 \
							$(OBJDIR_GRID)/Mod_Package_FileIO.$(OBJ) \
							$(OBJDIR_GRID)/Mod_Package_NACA.$(OBJ) \
							$(OBJDIR_GRID)/Mod_Package_Grid.$(OBJ) \
							$(OBJDIR_GRID)/Mod_Package_Flow.$(OBJ) \
							$(OBJDIR_GRID)/Mod_Package_ViewFlow.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<
$(OBJDIR_GRID)/Flow_Field_NACA.$(OBJ)		:	$(SRCDIR_GRID)/Flow_Field_NACA.f90 \
							$(OBJDIR_GRID)/Mod_Package_FileIO.$(OBJ) \
							$(OBJDIR_GRID)/Mod_Package_NACA.$(OBJ) \
							$(OBJDIR_GRID)/Mod_Package_Flow.$(OBJ) \
							$(OBJDIR_GRID)/Mod_Package_Icing.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<
$(OBJDIR_GRID)/View_Flow_NACA.$(OBJ)		:	$(SRCDIR_GRID)/View_Flow_NACA.f90 \
							$(OBJDIR_GRID)/Mod_Package_FileIO.$(OBJ) \
							$(OBJDIR_GRID)/Mod_Package_NACA.$(OBJ) \
							$(OBJDIR_GRID)/Mod_Package_Flow.$(OBJ) \
							$(OBJDIR_GRID)/Mod_Package_ViewFlow.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<

$(OBJDIR_MPSINI)/airfoil_2D.$(OBJ)		:	$(SRCDIR_MPSINI)/airfoil_2D.f90 \
							$(OBJDIR_MPSINI)/mod_is_2D.$(OBJ) \
							$(OBJDIR_MPSINI)/mod_is.$(OBJ) \
							$(OBJDIR_MPSINI)/mod_is_airfoil_2D.$(OBJ)

	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<

$(OBJDIR_MPS)/mps_icing_2D.$(OBJ)		:	$(SRCDIR_MPS)/mps_icing_2D.f90 \
							$(OBJDIR_MPS)/mod_mps_2D.$(OBJ) \
							$(OBJDIR_MPS)/mod_mps.$(OBJ) \
							$(OBJDIR_MPS)/mod_mps_icing_2D.$(OBJ) \
							$(OBJDIR_GRID)/Mod_Package_FileIO.$(OBJ) \
							$(OBJDIR_GRID)/Mod_Package_NACA.$(OBJ) \
							$(OBJDIR_GRID)/Mod_Package_Flow.$(OBJ) \
							$(OBJDIR_GRID)/Mod_Package_Grid.$(OBJ) \
							$(OBJDIR_MPS)/mod_package_droplet.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<

$(OBJDIR_ICM)/DeleteInnerParticle.$(OBJ)	:	$(SRCDIR_ICM)/DeleteInnerParticle.f90
	$(FC) $(FLAGS2) $(OPT_O) $@ $(OPT_C) $<

# モジュール
$(OBJDIR_GRID)/Mod_Package_FileIO.$(OBJ)		: $(MODDIR_GRID)/Mod_Package_FileIO.f90
	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<
$(OBJDIR_GRID)/Mod_Package_NACA.$(OBJ)			: $(MODDIR_GRID)/Mod_Package_NACA.f90
	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<
$(OBJDIR_GRID)/Mod_Package_Equation.$(OBJ)		: $(MODDIR_GRID)/Mod_Package_Equation.f90
	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<
$(OBJDIR_GRID)/Mod_Package_Grid.$(OBJ)			: $(MODDIR_GRID)/Mod_Package_Grid.f90
	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<
$(OBJDIR_GRID)/Mod_Package_Flow.$(OBJ)			: $(MODDIR_GRID)/Mod_Package_Flow.f90
	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<
$(OBJDIR_GRID)/Mod_Package_ViewFlow.$(OBJ)		: $(MODDIR_GRID)/Mod_Package_ViewFlow.f90
	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<
$(OBJDIR_GRID)/Mod_Package_Icing.$(OBJ)			: $(MODDIR_GRID)/Mod_Package_Icing.f90
	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<

$(OBJDIR_MPSINI)/mod_is.$(OBJ)				: $(MODDIR_MPSINI)/mod_is.f90
	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<
$(OBJDIR_MPSINI)/mod_is_2D.$(OBJ)			: $(MODDIR_MPSINI)/mod_is_2D.f90
	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<
$(OBJDIR_MPSINI)/mod_is_airfoil_2D.$(OBJ)		: $(MODDIR_MPSINI)/mod_is_airfoil_2D.f90
	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<

$(OBJDIR_MPS)/mod_mps.$(OBJ)				: $(MODDIR_MPS)/mod_mps.f90
	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<
$(OBJDIR_MPS)/mod_mps_2D.$(OBJ)				: $(MODDIR_MPS)/mod_mps_2D.f90
	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<
$(OBJDIR_MPS)/mod_mps_icing_2D.$(OBJ)			: $(MODDIR_MPS)/mod_mps_icing_2D.f90
	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<
$(OBJDIR_MPS)/mod_package_droplet.$(OBJ)		: $(MODDIR_MPS)/mod_package_droplet.f90
	$(FC) $(FLAGS) $(OPT_O) $@ $(OPT_C) $<

