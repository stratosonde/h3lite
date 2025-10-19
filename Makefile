# H3Lite Makefile
# Build system for h3lite library and tests

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
H3LITE_SRCS = $(SRC_DIR)/h3lite.c $(SRC_DIR)/h3lite_faceijk.c \
              $(SRC_DIR)/h3lite_basecells.c $(SRC_DIR)/h3lite_neighbor.c \
              $(SRC_DIR)/h3lite_regions_table.c

# For STM32 build (cross-compilation)
STM32_CC = arm-none-eabi-gcc
STM32_AR = arm-none-eabi-ar
STM32_CFLAGS = -mcpu=cortex-m4 -mthumb -mfloat-abi=hard -mfpu=fpv4-sp-d16 \
               -DUSE_HAL_DRIVER -DSTM32F411xE -Os -Wall -Wextra -I./include \
               -ffunction-sections -fdata-sections -fno-common
STM32_LDFLAGS = -lm --specs=nano.specs

# Create directories
$(shell mkdir -p $(BIN_DIR) $(OBJ_DIR))

# Default target
all: lib

# Static library build
lib: $(BIN_DIR)/libh3lite.a

$(BIN_DIR)/libh3lite.a: $(H3LITE_SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
	$(AR) rcs $@ $^

# Object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

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

# Nearest neighbor test
nearest_test: $(BIN_DIR)/h3lite_nearest_test

$(BIN_DIR)/h3lite_nearest_test: $(OBJ_DIR)/h3lite_nearest_test.o $(H3LITE_SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/h3lite_nearest_test.o: test/h3lite_nearest_test.c
	$(CC) $(CFLAGS) -c $< -o $@

# Run all tests
test: grid_test nearest_test
	@echo "\n=== Running Grid Coverage Test ==="
	./$(BIN_DIR)/h3lite_grid_test
	@echo "\n=== Running Nearest Neighbor Test ==="
	./$(BIN_DIR)/h3lite_nearest_test

# Generate lookup tables from GeoJSON region files
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


# Clean
clean:
	rm -rf $(BIN_DIR)/* $(OBJ_DIR)/* *.su

.PHONY: all lib stm32lib grid_test nearest_test test lookup_table clean stm32-size
