########################################################################################
######## マクロ定義								########
########################################################################################
FC = ifort
#FC = gfortran
OPT_O = -o
OPT_C = -c
#FLAGS = -fbounds-check -ffpe-trap=invalid,zero,overflow -fbacktrace
#FLAGS = -O2
#FLAGS = -fast
#FLAGS = -CB -fpe0 -traceback -g
#FLAGS = -mcmodel=large -shared-intel
FLAGS = -O2 -qopenmp -mcmodel=large
# -traceback
# -check bounds
#FLAGS = -fast -openmp
EXE = exe
OBJ = o
#MOD = ../Library/
DEL = rm
#DEL = del
ECHO   = echo
ALLEXE = Grid.$(EXE) \
         Overset.$(EXE) \
         IniFlow.$(EXE) \
         Flow.$(EXE) \
         ViewFlow.$(EXE) \
         Droplet.$(EXE) \
         Icing.$(EXE) \
         Validation.$(EXE) \
         Remesh.$(EXE)

########################################################################################
########## サフィックスルール							########
########################################################################################
.SUFFIXES: $(OBJ) .f90
.f90.o:
	$(FC) $(FLAGS) $(OPT_C) $< 

########################################################################################
########## 依存関係とコマンド							########
########################################################################################
# デフォルトターゲット
.PHONY: all
all: $(ALLEXE)

# 整理用
.PHONY: clean
clean : 
	$(DEL) *.$(OBJ) *.mod

# 実行ファイル
Grid.$(EXE)		: Grid_Generation_NACA.$(OBJ) \
	              	  Mod_Package_FileIO.$(OBJ) \
                  	  Mod_Package_NACA.$(OBJ) \
                  	  Mod_Package_Equation.$(OBJ) \
                 	  Mod_Package_Grid.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^
Overset.$(EXE)		: Overset_Interpolation_NACA.$(OBJ) \
                  	  Mod_Package_FileIO.$(OBJ) \
                  	  Mod_Package_NACA.$(OBJ) \
                  	  Mod_Package_Grid.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^
IniFlow.$(EXE)		: Initial_Flow_NACA.$(OBJ) \
                  	  Mod_Package_FileIO.$(OBJ) \
                  	  Mod_Package_NACA.$(OBJ) \
                  	  Mod_Package_Grid.$(OBJ) \
                  	  Mod_Package_Flow.$(OBJ) \
                  	  Mod_Package_ViewFlow.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^
Flow.$(EXE)		: Flow_Field_NACA.$(OBJ) \
                  	  Mod_Package_FileIO.$(OBJ) \
                  	  Mod_Package_NACA.$(OBJ) \
                  	  Mod_Package_Flow.$(OBJ) \
                 	  Mod_Package_Icing.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^
ViewFlow.$(EXE)		: View_Flow_NACA.$(OBJ) \
                  	  Mod_Package_FileIO.$(OBJ) \
                  	  Mod_Package_NACA.$(OBJ) \
                  	  Mod_Package_Flow.$(OBJ) \
                  	  Mod_Package_ViewFlow.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^
Droplet.$(EXE)		: Droplet_Trajectory_NACA.$(OBJ) \
                  	  Mod_Package_FileIO.$(OBJ) \
                  	  Mod_Package_NACA.$(OBJ) \
                  	  Mod_Package_Droplet.$(OBJ) \
                  	  Mod_Package_Grid.$(OBJ) \
                  	  Mod_Package_Flow.$(OBJ) \
                     	  mt95.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^
Icing.$(EXE)		: Ice_Accretion_NACA.$(OBJ) \
                  	  Mod_Package_FileIO.$(OBJ) \
                  	  Mod_Package_NACA.$(OBJ) \
                  	  Mod_Package_Flow.$(OBJ) \
                 	  Mod_Package_Icing.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^
Validation.$(EXE)	: Validation_NACA.$(OBJ) \
                    	  Mod_Package_FileIO.$(OBJ) \
                    	  Mod_Package_NACA.$(OBJ) \
                    	  Mod_Package_Flow.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^
Remesh.$(EXE)		: Remesh_NACA.$(OBJ) \
                  	  Mod_Package_FileIO.$(OBJ) \
                  	  Mod_Package_NACA.$(OBJ) \
                  	  Mod_Package_Equation.$(OBJ) \
                  	  Mod_Package_Grid.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^

# オブジェクト
Grid_Generation_NACA.$(OBJ)		: Grid_Generation_NACA.f90 \
                            		  Mod_Package_FileIO.$(OBJ) \
                            		  Mod_Package_NACA.$(OBJ) \
                            		  Mod_Package_Equation.$(OBJ) \
                                          Mod_Package_Grid.$(OBJ)
Overset_Interpolation_NACA.$(OBJ)	: Overset_Interpolation_NACA.f90 \
                            		  Mod_Package_FileIO.$(OBJ) \
                            		  Mod_Package_NACA.$(OBJ) \
                            		  Mod_Package_Grid.$(OBJ)
Initial_Flow_NACA.$(OBJ)		: Initial_Flow_NACA.f90 \
                            		  Mod_Package_FileIO.$(OBJ) \
                            		  Mod_Package_NACA.$(OBJ) \
                            		  Mod_Package_Grid.$(OBJ) \
                            		  Mod_Package_Flow.$(OBJ) \
                            		  Mod_Package_ViewFlow.$(OBJ)
Flow_Field_NACA.$(OBJ)			: Flow_Field_NACA.f90 \
                            		  Mod_Package_FileIO.$(OBJ) \
                            		  Mod_Package_NACA.$(OBJ) \
                            		  Mod_Package_Flow.$(OBJ) \
                            		  Mod_Package_Icing.$(OBJ)
View_Flow_NACA.$(OBJ)			: View_Flow_NACA.f90 \
                            		  Mod_Package_FileIO.$(OBJ) \
                            		  Mod_Package_NACA.$(OBJ) \
                            		  Mod_Package_Flow.$(OBJ) \
                            		  Mod_Package_ViewFlow.$(OBJ)
Droplet_Trajectory_NACA.$(OBJ)		: Droplet_Trajectory_NACA.f90 \
                            		  Mod_Package_FileIO.$(OBJ) \
                            		  Mod_Package_NACA.$(OBJ) \
                            		  Mod_Package_Droplet.$(OBJ) \
                            		  Mod_Package_Grid.$(OBJ) \
                            		  Mod_Package_Flow.$(OBJ) \
                                          mt95.$(OBJ)
Ice_Accretion_NACA.$(OBJ)		: Ice_Accretion_NACA.f90 \
                        		  Mod_Package_FileIO.$(OBJ) \
                           		  Mod_Package_NACA.$(OBJ) \
                           		  Mod_Package_Flow.$(OBJ) \
                            		  Mod_Package_Icing.$(OBJ)
Remesh_NACA.$(OBJ)			: Remesh_NACA.f90 \
                            		  Mod_Package_FileIO.$(OBJ) \
                            		  Mod_Package_NACA.$(OBJ) \
                            		  Mod_Package_Equation.$(OBJ) \
                                          Mod_Package_Grid.$(OBJ)
Validation_NACA.$(OBJ)			: Validation_NACA.f90 \
                            		  Mod_Package_FileIO.$(OBJ) \
                            		  Mod_Package_NACA.$(OBJ) \
                            		  Mod_Package_Flow.$(OBJ)


# モジュール
Mod_Package_FileIO.$(OBJ) 	: Mod_Package_FileIO.f90
Mod_Package_NACA.$(OBJ)   	: Mod_Package_NACA.f90
Mod_Package_Equation.$(OBJ)	: Mod_Package_Equation.f90
Mod_Package_Grid.$(OBJ)   	: Mod_Package_Grid.f90
Mod_Package_Flow.$(OBJ)   	: Mod_Package_Flow.f90
Mod_Package_ViewFlow.$(OBJ)   	: Mod_Package_ViewFlow.f90
Mod_Package_Droplet.$(OBJ)   	: Mod_Package_Droplet.f90
mt95.$(OBJ)  		 	: mt95.f90
Mod_Package_Icing.$(OBJ)  	: Mod_Package_Icing.f90
mod_icing.$(OBJ)  		: mod_icing.f90
mod_icing2.$(OBJ)  		: mod_icing2.f90

## オブジェクト
#Grid_Generation_NACA.$(OBJ)		: Grid_Generation_NACA.f90 \
#                            		$(MOD)Mod_Package_FileIO.$(OBJ) \
#                            		$(MOD)Mod_Package_NACA.$(OBJ) \
#                            		$(MOD)Mod_Package_Equation.$(OBJ) \
#                                       $(MOD)Mod_Package_Grid.$(OBJ)
#Overset_Interpolation_NACA.$(OBJ)	: Overset_Interpolation_NACA.f90 \
#                            		$(MOD)Mod_Package_FileIO.$(OBJ) \
#                            		$(MOD)Mod_Package_NACA.$(OBJ) \
#                            		$(MOD)Mod_Package_Grid.$(OBJ)
#Initial_Flow_NACA.$(OBJ)		: Initial_Flow_NACA.f90 \
#                            		$(MOD)Mod_Package_FileIO.$(OBJ) \
#                            		$(MOD)Mod_Package_NACA.$(OBJ) \
#                            		$(MOD)Mod_Package_Flow.$(OBJ)
#Flow_Field_NACA.$(OBJ)			: Flow_Field_NACA.f90 \
#                            		$(MOD)Mod_Package_FileIO.$(OBJ) \
#                            		$(MOD)Mod_Package_NACA.$(OBJ) \
#                            		$(MOD)Mod_Package_Flow.$(OBJ)
#View_Flow_NACA.$(OBJ)			: View_Flow_NACA.f90 \
#                            		$(MOD)Mod_Package_FileIO.$(OBJ) \
#                            		$(MOD)Mod_Package_NACA.$(OBJ) \
#                            		$(MOD)Mod_Package_ViewFlow.$(OBJ)
#Droplet_Trajectory_NACA.$(OBJ)		: Droplet_Trajectory_NACA.f90 \
#                            		$(MOD)Mod_Package_FileIO.$(OBJ) \
#                            		$(MOD)Mod_Package_NACA.$(OBJ) \
#                            		$(MOD)Mod_Package_Droplet.$(OBJ) \
#                                       $(MOD)mt95.$(OBJ)
#Ice_Accretion_NACA.$(OBJ)		: Ice_Accretion_NACA.f90 \
#                        		$(MOD)Mod_Package_FileIO.$(OBJ) \
#                           		$(MOD)Mod_Package_NACA.$(OBJ) \
#                            		$(MOD)Mod_Package_Icing.$(OBJ)
#Validation_NACA.$(OBJ)			: Validation_NACA.f90 \
#                        		$(MOD)Mod_Package_FileIO.$(OBJ) \
#                           		$(MOD)Mod_Package_NACA.$(OBJ) \
#                            		$(MOD)Mod_Package_Flow.$(OBJ)
#Remesh_NACA.$(OBJ)			: Remesh_NACA.f90 \
#                            		$(MOD)Mod_Package_FileIO.$(OBJ) \
#                            		$(MOD)Mod_Package_NACA.$(OBJ) \
#                            		$(MOD)Mod_Package_Equation.$(OBJ) \
#                                       $(MOD)Mod_Package_Grid.$(OBJ)
#
## モジュール
#$(MOD)Mod_Package_FileIO.$(OBJ) 	: $(MOD)Mod_Package_FileIO.f90
#$(MOD)Mod_Package_NACA.$(OBJ)   	: $(MOD)Mod_Package_NACA.f90
#$(MOD)Mod_Package_Equation.$(OBJ)	: $(MOD)Mod_Package_Equation.f90
#$(MOD)Mod_Package_Grid.$(OBJ)   	: $(MOD)Mod_Package_Grid.f90
#$(MOD)Mod_Package_Flow.$(OBJ)   	: $(MOD)Mod_Package_Flow.f90
#$(MOD)Mod_Package_ViewFlow.$(OBJ)   	: $(MOD)Mod_Package_ViewFlow.f90
#$(MOD)Mod_Package_Droplet.$(OBJ)   	: $(MOD)Mod_Package_Droplet.f90
#$(MOD)mt95.$(OBJ)  		 	: $(MOD)mt95.f90
#$(MOD)Mod_Package_Icing.$(OBJ)  	: $(MOD)Mod_Package_Icing.f90
#$(MOD)Mod_Package_Validation.$(OBJ)  	: $(MOD)Mod_Package_Validation.f90
