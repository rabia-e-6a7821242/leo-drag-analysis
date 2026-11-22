FC      = gfortran
FFLAGS  = -O2 -Wall -std=f2008 -fimplicit-none
SRC     = src
BUILD   = build
BIN     = bin

MODS    = constants atmosphere space_weather spacecraft
MOD_OBJ = $(addprefix $(BUILD)/,$(addsuffix .o,$(MODS)))

all: $(BIN)/orbit_decay $(BIN)/density_profile

$(BIN)/%: $(BUILD)/%.o $(MOD_OBJ) | $(BIN)
	$(FC) $(FFLAGS) -o $@ $(MOD_OBJ) $<

$(BUILD)/%.o: $(SRC)/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -J$(BUILD) -c $< -o $@

# module dependencies
$(BUILD)/atmosphere.o:      $(BUILD)/constants.o
$(BUILD)/space_weather.o:   $(BUILD)/constants.o
$(BUILD)/spacecraft.o:      $(BUILD)/constants.o
$(BUILD)/orbit_decay.o:     $(MOD_OBJ)
$(BUILD)/density_profile.o: $(BUILD)/constants.o $(BUILD)/atmosphere.o

$(BUILD) $(BIN):
	mkdir -p $@

clean:
	rm -rf $(BUILD) $(BIN)

.PHONY: all clean
.PRECIOUS: $(BUILD)/%.o
