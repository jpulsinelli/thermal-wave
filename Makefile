MACHINE = juliana

OPT_LEVEL = OPTIMIZE
FLAGS     = $(FLAGS_$(OPT_LEVEL))

include ./Makefile_Main

all: Thermal

Thermal: \
	$(objects) ThermalMain.o
	$(FLINKER) $(FLAGS) -o Thermal_$(MACHINE) \
	$(objects) ThermalMain.o $(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f Thermal_$(MACHINE)
