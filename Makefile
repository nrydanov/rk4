# Переменные для путей
SIM_DIR = simulation
BUILD_DIR = $(SIM_DIR)/build
DISPLAY_DIR = display
CONFIG = config.yaml
RESULTS = results.csv

.PHONY: all run build plot clean

all: execute plot

release:
	@echo "--- Building C++ simulation (Release) ---"
	rm -rf $(BUILD_DIR)
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && cmake -DCMAKE_BUILD_TYPE=Release .. && cmake --build . -j$(nproc)

build:
	@echo "--- Building C++ simulation ---"
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && cmake ..
	cmake --build $(BUILD_DIR)

execute: build
	@echo "--- Running simulation ---"
	./$(BUILD_DIR)/vdp_sim $(SIM_DIR)/$(CONFIG) -o ./$(RESULTS)

plot:
	@echo "--- Displaying results ---"
	cd $(DISPLAY_DIR) && uv run display.py ../$(RESULTS)

clean:
	rm -f $(RESULTS)
	rm -f *.png
	rm -rf $(BUILD_DIR)
	@echo "Cleaned up results and plots."
