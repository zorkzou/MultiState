#
# This Makefile may be used for gfortran under Linux or Windows (msys or cygwin).
#
# Written by W. Zou <qcband@gmail.com> on Oct. 31, 2022.
#

DIR        = ../
NAME       = multistate.exe
EXE        = $(DIR)$(NAME)

OBJ        = object

FCOM       = gfortran
COMPFLAGS  = -O3
#            Looking for .o & .mod files here!
COMPFLAGS += -I$(OBJ) -J$(OBJ)

SRC90      = multistate.f90
SRC77      = blas.f  lapack.f
MODS       =

F90OBJ     = $(addprefix $(OBJ)/, $(SRC90:.f90=.o))
F77OBJ     = $(addprefix $(OBJ)/, $(SRC77:.f=.o))
MODOBJ     = $(addprefix $(OBJ)/, $(MODS:.f90=.o))

all: $(MODOBJ) $(F90OBJ) $(F77OBJ)
	@echo
	$(FCOM) $(OBJ)/*.o $(COMPFLAGS) -o $(EXE)
	@echo
	@echo "  <<< $(NAME) has been compiled successfully! >>>"
	@echo

$(MODOBJ): $(OBJ)/%.o: %.f90
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f90

$(F90OBJ): $(OBJ)/%.o: %.f90
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f90

$(F77OBJ): $(OBJ)/%.o: %.f
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f

.PHONY: clean
clean:
	rm -f $(OBJ)/*.o $(OBJ)/*.mod $(OBJ)/*.MOD
