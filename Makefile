# H3Lite Makefile for testing and building the library
# Optimized for minimal code size

CC = gcc
AR = ar
CFLAGS = -Wall -Wextra -O3 -g -I./include
LDFLAGS = -lm

# Directories
SRC_DIR = src
INC_DIR = include
BIN_DIR = bin
OBJ_DIR = obj

# Source files
H3LITE_SRCS = $(SRC_DIR)/h3lite.c $(SRC_DIR)/h3lite_faceijk.c $(SRC_DIR)/h3lite_basecells.c $(SRC_DIR)/h3lite_neighbor.c $(SRC_DIR)/h3lite_regions_table.c

# For STM32 build (cross-compilation)
STM32_CC = arm-none-eabi-gcc
STM32_AR = arm-none-eabi-ar
STM32_CFLAGS = -mcpu=cortex-m4 -mthumb -mfloat-abi=hard -mfpu=fpv4-sp-d16 \
               -DUSE_HAL_DRIVER -DSTM32F411xE -Os -Wall -Wextra -I./include \
               -ffunction-sections -fdata-sections -fno-common
STM32_LDFLAGS = -lm --specs=nano.specs

# Create directories
$(shell mkdir -p $(BIN_DIR) $(OBJ_DIR))

# Default build
all: grid_test

# Static library build
lib: $(BIN_DIR)/libh3lite.a

$(BIN_DIR)/libh3lite.a: $(H3LITE_SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
	$(AR) rcs $@ $^

# STM32 static library build
stm32lib: $(H3LITE_SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%_stm32.o)
	$(STM32_AR) rcs $(BIN_DIR)/libh3lite_stm32.a $^

$(OBJ_DIR)/%_stm32.o: $(SRC_DIR)/%.c
	$(STM32_CC) $(STM32_CFLAGS) -c $< -o $@

# Test programs
grid_test: $(BIN_DIR)/h3lite_grid_test

$(BIN_DIR)/h3lite_grid_test: $(OBJ_DIR)/h3lite_grid_test.o $(H3LITE_SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/h3lite_grid_test.o: test/h3lite_grid_test.c
	$(CC) $(CFLAGS) -c $< -o $@

# Nearest regions test
nearest_test: $(BIN_DIR)/h3lite_nearest_test

$(BIN_DIR)/h3lite_nearest_test: $(OBJ_DIR)/h3lite_nearest_test.o $(H3LITE_SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/h3lite_nearest_test.o: test/h3lite_nearest_test.c
	$(CC) $(CFLAGS) -c $< -o $@

# Region lookup example
region_example: $(BIN_DIR)/region_lookup_example

$(BIN_DIR)/region_lookup_example: $(OBJ_DIR)/region_lookup_example.o $(H3LITE_SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/region_lookup_example.o: test/region_lookup_example.c
	$(CC) $(CFLAGS) -c $< -o $@

# Debug San Francisco test
debug_sf: $(BIN_DIR)/debug_sf

$(BIN_DIR)/debug_sf: $(OBJ_DIR)/debug_sf.o $(H3LITE_SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/debug_sf.o: test/debug_sf.c
	$(CC) $(CFLAGS) -c $< -o $@

# Debug overage handling
debug_overage: $(BIN_DIR)/debug_overage

$(BIN_DIR)/debug_overage: $(OBJ_DIR)/debug_overage.o $(H3LITE_SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/debug_overage.o: test/debug_overage.c
	$(CC) $(CFLAGS) -c $< -o $@

# H3 comparison tests (two-stage approach to avoid linking conflicts)
compare: $(BIN_DIR)/h3_comparison_test $(BIN_DIR)/h3lite_compare

compare-run: compare
	@echo "Running H3 comparison tests..."
	@cd $(BIN_DIR) && ./h3_comparison_test && ./h3lite_compare

$(BIN_DIR)/h3_comparison_test: $(OBJ_DIR)/h3_comparison_test.o
	$(CC) $(CFLAGS) -o $@ $^ -I/home/englotk/hplans/h3/build/src/h3lib/include -L/home/englotk/hplans/h3/build/lib -lh3 -Wl,-rpath=/home/englotk/hplans/h3/build/lib $(LDFLAGS)

# Object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# Run all tests
test_all: grid_test nearest_test
	@echo "\n=== Running Grid Coverage Test ==="
	./$(BIN_DIR)/h3lite_grid_test
	@echo "\n=== Running Nearest Regions Test ==="
	./$(BIN_DIR)/h3lite_nearest_test

$(OBJ_DIR)/h3_comparison_test.o: test/h3_comparison_test.c
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN_DIR)/h3lite_compare: $(OBJ_DIR)/h3lite_compare.o $(H3LITE_SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/h3lite_compare.o: test/h3lite_compare.c
	$(CC) $(CFLAGS) -c $< -o $@

# Generate lookup table
lookup_table:
	python3 generate_lookup_table.py

# Memory usage analysis for STM32
stm32-size: stm32lib
	@echo "=== Flash and RAM Usage Analysis ==="
	@echo "Total library size:"
	@arm-none-eabi-size $(BIN_DIR)/libh3lite_stm32.a
	@echo "\nSize by object file:"
	@arm-none-eabi-size $(OBJ_DIR)/*_stm32.o
	@echo "\nDetailed section breakdown:"
	@arm-none-eabi-size -A $(BIN_DIR)/libh3lite_stm32.a
	@echo "\nLargest functions (top 20):"
	@arm-none-eabi-nm --print-size --size-sort $(BIN_DIR)/libh3lite_stm32.a | grep -v ' [a-z] ' | tail -20
	@echo "\nRAM usage (data + BSS sections):"
	@arm-none-eabi-nm --print-size --size-sort $(BIN_DIR)/libh3lite_stm32.a | grep -E " [bdvBDV] " | sort -k2nr

# Test executable for memory usage analysis
stm32-test: $(OBJ_DIR)/stm32_size_test_stm32.o $(H3LITE_SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%_stm32.o)
	$(STM32_CC) $(STM32_CFLAGS) -o $(BIN_DIR)/stm32_test.elf $^ $(STM32_LDFLAGS) -Wl,-Map=$(BIN_DIR)/stm32_test.map
	@echo "\n=== Complete STM32 Executable Size ==="
	@arm-none-eabi-size $(BIN_DIR)/stm32_test.elf
	@echo "\nSection Headers:"
	@arm-none-eabi-objdump -h $(BIN_DIR)/stm32_test.elf

$(OBJ_DIR)/stm32_size_test_stm32.o: test/stm32_size_test.c
	$(STM32_CC) $(STM32_CFLAGS) -c $< -o $@

# Stack usage analysis (requires -fstack-usage compiler flag)
stm32-stack: stm32lib
	@echo "=== Stack Usage Analysis ==="
	@echo "Adding -fstack-usage to compilation flags for analysis"
	@$(STM32_CC) $(STM32_CFLAGS) -fstack-usage -c $(SRC_DIR)/h3lite.c -o $(OBJ_DIR)/h3lite_stack.o
	@$(STM32_CC) $(STM32_CFLAGS) -fstack-usage -c $(SRC_DIR)/h3lite_faceijk.c -o $(OBJ_DIR)/h3lite_faceijk_stack.o
	@echo "\nStack usage by function:"
	@cat *.su 2>/dev/null || echo "No stack usage files found"

# Clean
clean:
	rm -rf $(BIN_DIR)/* $(OBJ_DIR)/* *.su

.PHONY: all lib stm32lib test compare lookup_table clean stm32-size stm32-test stm32-stack
