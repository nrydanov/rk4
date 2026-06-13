# Переменные для путей
SIM_DIR = simulation
BUILD_DIR = $(SIM_DIR)/build
DISPLAY_DIR = display
RESULTS     = results.csv
RESULTS_PQ  = results.parquet
CONFIG = config.yaml

.PHONY: demo run build plot clean sync_energy phase_slopes multistability validation to-parquet

# Быстрая демонстрация: сборка + базовая симуляция + просмотр траекторий и фазовых портретов
demo: release run-demo plot-demo

# То же, что demo, но без открытия окна — только сохраняет картинки
save: release run-demo save-demo

# Эксперимент: карта когерентности на плоскости (δ1, δ2) для набора ε
sync_energy: release run-sync-energy to-parquet plot-heatmap

# Эксперимент: карта когерентности + режимы фазовой синхронизации на плоскости (δ1, δ2)
phase_slopes: release run-phase-slopes to-parquet plot-phase-slopes

# Эксперимент: мультистабильность — 4 варианта связи, случайные начальные условия
multistability: release run-multistability to-parquet

# Валидация пайплайна против аналитики: амплитуда/частота ВдП, порядок RK4,
# полная синхронизация идентичной пары, скейлинг языка Адлера, симметрия звезды
validation: release
	@echo "--- Running validation suite ---"
	./$(BUILD_DIR)/vdp_validation

# Запуск базовой симуляции (одна траектория)
run-demo:
	@echo "--- Running simulation ---"
	./$(BUILD_DIR)/vdp_sim $(SIM_DIR)/$(CONFIG) -o ./$(RESULTS)

# Запуск эксперимента по когерентности: перебор (δ1, δ2, ε), результат — results.csv
run-sync-energy:
	@echo "--- Running sync energy experiment ---"
	./$(BUILD_DIR)/vdp_sync_energy $(SIM_DIR)/experiments/sync_energy/$(CONFIG) -o ./$(RESULTS)

# Запуск эксперимента по фазовым наклонам: перебор (δ1, δ2, ε), результат — results.csv
run-phase-slopes:
	@echo "--- Running phase slopes experiment ---"
	./$(BUILD_DIR)/vdp_phase_slopes $(SIM_DIR)/experiments/phase_slopes/$(CONFIG) -o ./$(RESULTS)

# Запуск эксперимента по мультистабильности: 4 варианта связи × случайные НУ
# Привязка потоков к ядрам (close/cores) — стабильная кэш-локальность на x86/Linux
run-multistability:
	@echo "--- Running multistability experiment ---"
	OMP_PROC_BIND=close OMP_PLACES=cores ./$(BUILD_DIR)/vdp_multistability $(SIM_DIR)/experiments/multistability/$(CONFIG) -o ./$(RESULTS)

# Конвертация results.csv → results.parquet (zstd, без потери точности)
to-parquet:
	@echo "--- Converting $(RESULTS) → $(RESULTS_PQ) ---"
	cd $(DISPLAY_DIR) && uv run python -c "\
import pandas as pd; \
df = pd.read_csv('../$(RESULTS)', comment='#'); \
df.to_parquet('../$(RESULTS_PQ)', compression='zstd', index=False); \
print(f'  {len(df):,} rows → $(RESULTS_PQ)')"

# Отрисовка траекторий и фазовых портретов с открытием окна
plot-demo:
	@echo "--- Displaying results ---"
	cd $(DISPLAY_DIR) && uv run display.py ../$(RESULTS) --show

# То же, что plot-demo, без открытия окна
save-demo:
	@echo "--- Saving results ---"
	cd $(DISPLAY_DIR) && uv run display.py ../$(RESULTS)

# Тепловая карта когерентности L(δ1, δ2) для каждого ε
plot-heatmap:
	@echo "--- Displaying heatmap results ---"
	cd $(DISPLAY_DIR) && uv run heatmap.py ../$(RESULTS)

# Парные картинки: когерентность + режимы синхронизации для каждого ε
plot-phase-slopes:
	@echo "--- Displaying phase sync results ---"
	cd $(DISPLAY_DIR) && uv run phase_sync.py ../$(RESULTS)

# Сборка всех бинарников в режиме Release с оптимизациями
release:
	@echo "--- Building C++ simulation (Release) ---"
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && \
		cmake -DCMAKE_BUILD_TYPE=Release .. && \
		cmake --build . -j$(nproc)

# Сборка в режиме Debug — для отладки, без оптимизаций
debug:
	@echo "--- Building C++ simulation (Debug) ---"
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && \
		cmake .. && \
		cmake --build . -j$(nproc)

# Удаление результатов, картинок и артефактов сборки
clean:
	rm -f $(RESULTS) $(RESULTS_PQ)
	rm -f *.png
	rm -rf $(BUILD_DIR)
	@echo "Cleaned up results and plots."
