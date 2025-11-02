########################################################################################
######## gfortran専用設定 ##############################################################
########################################################################################
.DEFAULT_GOAL := all

FC      = gfortran
EXE     = exe
OBJ     = o
MODDIR  = mod

FLAGS   = -O2 -fopenmp -J$(MODDIR) -I$(MODDIR)
OPT_O   = -o
OPT_C   = -c
DEL     = rm -f

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
######## ビルドルール ##################################################################
########################################################################################
# .mod 出力先を先に用意
$(MODDIR):
	mkdir -p $(MODDIR)

# すべての .f90 → .o （order-only で mod/ を先に作る）
%.o: %.f90 | $(MODDIR)
	$(FC) $(FLAGS) $(OPT_C) $<

.PHONY: all clean
all: $(ALLEXE)

clean:
	$(DEL) *.$(OBJ) *.mod
	$(DEL) -r $(MODDIR) 2>/dev/null || true
	# $(DEL) $(ALLEXE)  # 実行ファイルも消したい場合は有効化

# ===== 実行ファイル（リンク） =====
Grid.$(EXE): Grid_Generation_NACA.$(OBJ) Mod_Package_FileIO.$(OBJ) Mod_Package_NACA.$(OBJ) Mod_Package_Equation.$(OBJ) Mod_Package_Grid.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^

Overset.$(EXE): Overset_Interpolation_NACA.$(OBJ) Mod_Package_FileIO.$(OBJ) Mod_Package_NACA.$(OBJ) Mod_Package_Grid.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^

IniFlow.$(EXE): Initial_Flow_NACA.$(OBJ) Mod_Package_FileIO.$(OBJ) Mod_Package_NACA.$(OBJ) Mod_Package_Grid.$(OBJ) Mod_Package_Flow.$(OBJ) Mod_Package_ViewFlow.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^

Flow.$(EXE): Flow_Field_NACA.$(OBJ) Mod_Package_FileIO.$(OBJ) Mod_Package_NACA.$(OBJ) Mod_Package_Flow.$(OBJ) Mod_Package_Icing.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^

ViewFlow.$(EXE): View_Flow_NACA.$(OBJ) Mod_Package_FileIO.$(OBJ) Mod_Package_NACA.$(OBJ) Mod_Package_Flow.$(OBJ) Mod_Package_ViewFlow.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^

Droplet.$(EXE): Droplet_Trajectory_NACA.$(OBJ) Mod_Package_FileIO.$(OBJ) Mod_Package_NACA.$(OBJ) Mod_Package_Droplet.$(OBJ) Mod_Package_Grid.$(OBJ) Mod_Package_Flow.$(OBJ) mt95.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^

Icing.$(EXE): Ice_Accretion_NACA.$(OBJ) Mod_Package_FileIO.$(OBJ) Mod_Package_NACA.$(OBJ) Mod_Package_Flow.$(OBJ) Mod_Package_Icing.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^

Validation.$(EXE): Validation_NACA.$(OBJ) Mod_Package_FileIO.$(OBJ) Mod_Package_NACA.$(OBJ) Mod_Package_Flow.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^

Remesh.$(EXE): Remesh_NACA.$(OBJ) Mod_Package_FileIO.$(OBJ) Mod_Package_NACA.$(OBJ) Mod_Package_Equation.$(OBJ) Mod_Package_Grid.$(OBJ)
	$(FC) $(FLAGS) $(OPT_O) $@ $^

# ===== 依存（ビルド順ヒント） =====
Grid_Generation_NACA.$(OBJ): Grid_Generation_NACA.f90 Mod_Package_FileIO.$(OBJ) Mod_Package_NACA.$(OBJ) Mod_Package_Equation.$(OBJ) Mod_Package_Grid.$(OBJ)
Overset_Interpolation_NACA.$(OBJ): Overset_Interpolation_NACA.f90 Mod_Package_FileIO.$(OBJ) Mod_Package_NACA.$(OBJ) Mod_Package_Grid.$(OBJ)
Initial_Flow_NACA.$(OBJ): Initial_Flow_NACA.f90 Mod_Package_FileIO.$(OBJ) Mod_Package_NACA.$(OBJ) Mod_Package_Grid.$(OBJ) Mod_Package_Flow.$(OBJ) Mod_Package_ViewFlow.$(OBJ)
Flow_Field_NACA.$(OBJ): Flow_Field_NACA.f90 Mod_Package_FileIO.$(OBJ) Mod_Package_NACA.$(OBJ) Mod_Package_Flow.$(OBJ) Mod_Package_Icing.$(OBJ)
View_Flow_NACA.$(OBJ): View_Flow_NACA.f90 Mod_Package_FileIO.$(OBJ) Mod_Package_NACA.$(OBJ) Mod_Package_Flow.$(OBJ) Mod_Package_ViewFlow.$(OBJ)
Droplet_Trajectory_NACA.$(OBJ): Droplet_Trajectory_NACA.f90 Mod_Package_FileIO.$(OBJ) Mod_Package_NACA.$(OBJ) Mod_Package_Droplet.$(OBJ) Mod_Package_Grid.$(OBJ) Mod_Package_Flow.$(OBJ) mt95.$(OBJ)
Ice_Accretion_NACA.$(OBJ): Ice_Accretion_NACA.f90 Mod_Package_FileIO.$(OBJ) Mod_Package_NACA.$(OBJ) Mod_Package_Flow.$(OBJ) Mod_Package_Icing.$(OBJ)
Remesh_NACA.$(OBJ): Remesh_NACA.f90 Mod_Package_FileIO.$(OBJ) Mod_Package_NACA.$(OBJ) Mod_Package_Equation.$(OBJ) Mod_Package_Grid.$(OBJ)
Validation_NACA.$(OBJ): Validation_NACA.f90 Mod_Package_FileIO.$(OBJ) Mod_Package_NACA.$(OBJ) Mod_Package_Flow.$(OBJ)
