# Переменные для путей
SIM_DIR = simulation
BUILD_DIR = $(SIM_DIR)/build
DISPLAY_DIR = display
RESULTS = results.csv

.PHONY: demo run build plot clean

demo: release run-demo plot-demo
heatmap: release run-heatmap plot-heatmap

plot-demo:
	@echo "--- Displaying results ---"
	cd $(DISPLAY_DIR) && uv run display.py ../$(RESULTS)

run-demo:
	@echo "--- Running simulation ---"
	./$(BUILD_DIR)/vdp_sim $(SIM_DIR)/$(CONFIG) -o ./$(RESULTS)

plot-heatmap:
	@echo "--- Displaying heatmap results ---"
	cd $(DISPLAY_DIR) && uv run heatmap.py ../$(RESULTS)

run-heatmap:
	@echo "--- Running heatmap generation ---"
	./$(BUILD_DIR)/vdp_heatmap $(SIM_DIR)/experiments/heatmap/$(CONFIG) -o ./$(RESULTS)

release:
	@echo "--- Building C++ simulation (Release) ---"
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && \
		CXX=g++-15 CC=gcc-15 \
		cmake -DCMAKE_BUILD_TYPE=Release .. && \
		cmake --build . -j$(nproc)

debug:
	@echo "--- Building C++ simulation (Debug) ---"
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && \
		CXX=g++-15 CC=gcc-15 \
		cmake .. && \
		cmake --build . -j$(nproc)

clean:
	rm -f $(RESULTS)
	rm -f *.png
	rm -rf $(BUILD_DIR)
	@echo "Cleaned up results and plots."
