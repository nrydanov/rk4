# Переменные для путей
SIM_DIR = simulation
BUILD_DIR = $(SIM_DIR)/build
DISPLAY_DIR = display
CONFIG = config.yaml
RESULTS = results.csv

.PHONY: all run build plot clean

all: execute plot

build:
	@echo "--- Building C++ simulation ---"
	cd $(SIM_DIR) && cmake --build ./build

execute: build
	@echo "--- Running simulation ---"
	./$(BUILD_DIR)/vdp_sim $(SIM_DIR)/$(CONFIG) -o ./$(RESULTS)

plot:
	@echo "--- Displaying results ---"
	cd $(DISPLAY_DIR) && uv run display.py ../$(RESULTS)

clean:
	rm -f $(RESULTS)
	rm -f *.png
	@echo "Cleaned up results and plots."
