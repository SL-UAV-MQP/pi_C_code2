# ============================================================================
# Unified Makefile for MQP Full System (Raspberry Pi 5B Deployment)
# ============================================================================
#
# Target: ARM Cortex-A76 (Raspberry Pi 5B)
# Dependencies: OpenBLAS, FFTW3, LibIIO (optional)
# C Standard: C99
#
# Build targets:
#   make all          - Build complete library, app, and tests
#   make app          - Build mqp_localize application
#   make lib          - Build static library
#   make test         - Build and run MUSIC unit tests
#   make test-sdr     - Build and run SDR tests (requires 3x PLUTO hardware)
#   make clean        - Remove build artifacts
#   make install      - Install to system (requires sudo)
#   make check-deps   - Verify dependencies installed
#   make info         - Display compiler and configuration info
#
# Options:
#   DEBUG=1           - Build with debug symbols, verbose, profiling
#   USE_SDR=1         - Include ADALM-PLUTO SDR support (LibIIO)
#
# ============================================================================

# Compiler and archiver
CC = gcc
AR = ar

# Base CFLAGS: C99, warnings, optimization for RPi 5B (ARMv8.2-A Cortex-A76)
CFLAGS = -std=c99 -Wall -Wextra -O3
CFLAGS += -march=armv8.2-a+fp16+dotprod -mtune=cortex-a76
CFLAGS += -ffast-math -fno-finite-math-only -funroll-loops -ftree-vectorize
CFLAGS += -I./include

# Dependency generation
DEPFLAGS = -MMD -MP

# Debug build (use: make DEBUG=1)
ifdef DEBUG
CFLAGS += -g -O0 -DMUSIC_VERBOSE=1 -DMUSIC_PROFILE=1
else
CFLAGS += -DNDEBUG
endif

# ARM NEON SIMD optimizations (enabled by default for RPi5)
CFLAGS += -DUSE_ARM_NEON=1

# ============================================================================
# Libraries
# ============================================================================

# Core: OpenBLAS (BLAS/LAPACK), FFTW3, math, pthreads
LDLIBS = -lopenblas -lfftw3 -lm -lpthread

# Optional: LibIIO for ADALM-PLUTO SDR
ifdef USE_SDR
CFLAGS += -DUSE_SDR=1
LDLIBS += -liio
SDR_SOURCES = $(SRC_DIR)/sdr_pluto.c
else
SDR_SOURCES =
endif

# ============================================================================
# Directories
# ============================================================================

SRC_DIR = src
INC_DIR = include
TEST_DIR = test
BUILD_DIR = build
LIB_DIR = lib

# ============================================================================
# Source Files Organization
# ============================================================================

# Core MUSIC AOA Pipeline + State Machine
MUSIC_SOURCES = $(SRC_DIR)/music_uca_6.c \
                $(SRC_DIR)/array_geometry.c \
                $(SRC_DIR)/steering_vector.c \
                $(SRC_DIR)/covariance.c \
                $(SRC_DIR)/eigendecomp.c \
                $(SRC_DIR)/spectrum.c \
                $(SRC_DIR)/peak_detection.c \
                $(SRC_DIR)/cfar_2d.c \
                $(SRC_DIR)/common_types.c \
                $(SRC_DIR)/state_machine.c

# Signal Processing
SIGNAL_SOURCES = $(SRC_DIR)/signal_processor.c \
                 $(SRC_DIR)/beamforming.c \
                 $(SRC_DIR)/noise_reduction.c

# Protocol Detection
PROTOCOL_SOURCES = $(SRC_DIR)/lte_detector.c \
                   $(SRC_DIR)/p25_detector.c \
                   $(SRC_DIR)/cell_tower_db.c

# Localization & Navigation
LOCALIZATION_SOURCES = $(SRC_DIR)/pathfinding.c \
                       $(SRC_DIR)/field_generation.c \
                       $(SRC_DIR)/map_fitting.c \
                       $(SRC_DIR)/aoa_triangulation.c

# Testing & Diagnostics
DIAG_SOURCES = $(SRC_DIR)/diagnostic_logger.c \
               $(SRC_DIR)/wpi_test_config.c

# All sources combined
ALL_SOURCES = $(MUSIC_SOURCES) $(SIGNAL_SOURCES) \
              $(PROTOCOL_SOURCES) $(LOCALIZATION_SOURCES) \
              $(DIAG_SOURCES) $(SDR_SOURCES)

# Object files
OBJECTS = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(ALL_SOURCES))

# Library output
LIB_NAME = libmqp_full.a
LIB_PATH = $(LIB_DIR)/$(LIB_NAME)

# Test executables
TEST_MUSIC_SRC = $(TEST_DIR)/test_music.c
TEST_MUSIC_BIN = $(BUILD_DIR)/test_music
TEST_SDR_SRC = $(TEST_DIR)/test_sdr_pluto.c
TEST_SDR_BIN = $(BUILD_DIR)/test_sdr_pluto
TEST_SIGPROC_SRC = $(TEST_DIR)/test_signal_processor.c
TEST_SIGPROC_BIN = $(BUILD_DIR)/test_signal_processor
TEST_AOA_SRC = $(TEST_DIR)/test_aoa_triangulation.c
TEST_AOA_BIN = $(BUILD_DIR)/test_aoa_triangulation
TEST_MODULES_SRC = $(TEST_DIR)/test_modules.c
TEST_MODULES_BIN = $(BUILD_DIR)/test_modules
TEST_WPI_SRC = $(TEST_DIR)/test_wpi_integration.c
TEST_WPI_BIN = $(BUILD_DIR)/test_wpi_integration

# Main application
MAIN_SRC = $(SRC_DIR)/main.c
MAIN_BIN = $(BUILD_DIR)/mqp_localize

# Installation paths
PREFIX = /usr/local
INSTALL_LIB_DIR = $(PREFIX)/lib
INSTALL_INC_DIR = $(PREFIX)/include/mqp

# ============================================================================
# Build Rules
# ============================================================================

.PHONY: all app lib test test-all test-sdr test-sigproc test-aoa test-modules test-wpi clean install uninstall help info check-deps run
.PHONY: lib-music lib-signal lib-protocol lib-localization

all: lib app test

# Create build directories
$(BUILD_DIR) $(LIB_DIR):
	mkdir -p $@

# Compile source files (with dependency generation)
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	@echo "  CC    $<"
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

# Include auto-generated dependency files
-include $(OBJECTS:.o=.d)

# Build static library
lib: $(LIB_PATH)

$(LIB_PATH): $(OBJECTS) | $(LIB_DIR)
	@echo "  AR    $(LIB_NAME) ($(words $(OBJECTS)) objects)"
	$(AR) rcs $@ $(OBJECTS)
	@echo "Library built: $@"

# ============================================================================
# Main Application
# ============================================================================

# Build main localization application
app: $(MAIN_BIN)

$(MAIN_BIN): $(MAIN_SRC) $(LIB_PATH) | $(BUILD_DIR)
	@echo "  CC    mqp_localize"
	$(CC) $(CFLAGS) $< -o $@ -L$(LIB_DIR) -lmqp_full $(LDLIBS)
	@echo "Application built: $@"

# Run application (convenience target)
run: $(MAIN_BIN)
	@echo ""
	@echo "Starting MQP Localization System..."
	@echo "======================================"
	$(MAIN_BIN) -s wpi -d

# ============================================================================
# Test Targets
# ============================================================================

# MUSIC unit tests (no SDR required)
test: $(TEST_MUSIC_BIN)
	@echo ""
	@echo "Running MUSIC unit tests..."
	@echo "======================================"
	$(TEST_MUSIC_BIN)

$(TEST_MUSIC_BIN): $(TEST_MUSIC_SRC) $(LIB_PATH) | $(BUILD_DIR)
	@echo "  CC    test_music"
	$(CC) $(CFLAGS) $< -o $@ -L$(LIB_DIR) -lmqp_full $(LDLIBS)

# SDR integration tests (requires USE_SDR=1 and 3x PLUTO hardware)
test-sdr: $(TEST_SDR_BIN)
	@echo ""
	@echo "Running SDR integration tests..."
	@echo "======================================"
	@echo "NOTE: Requires 3x ADALM-PLUTO SDR connected via USB"
	$(TEST_SDR_BIN)

$(TEST_SDR_BIN): $(TEST_SDR_SRC) $(LIB_PATH) | $(BUILD_DIR)
	@echo "  CC    test_sdr_pluto"
	$(CC) $(CFLAGS) $< -o $@ -L$(LIB_DIR) -lmqp_full $(LDLIBS)

# Signal processor tests
test-sigproc: $(TEST_SIGPROC_BIN)
	@echo ""
	@echo "Running signal processor tests..."
	@echo "======================================"
	$(TEST_SIGPROC_BIN)

$(TEST_SIGPROC_BIN): $(TEST_SIGPROC_SRC) $(LIB_PATH) | $(BUILD_DIR)
	@echo "  CC    test_signal_processor"
	$(CC) $(CFLAGS) $< -o $@ -L$(LIB_DIR) -lmqp_full $(LDLIBS)

# AOA triangulation tests
test-aoa: $(TEST_AOA_BIN)
	@echo ""
	@echo "Running AOA triangulation tests..."
	@echo "======================================"
	$(TEST_AOA_BIN)

$(TEST_AOA_BIN): $(TEST_AOA_SRC) $(LIB_PATH) | $(BUILD_DIR)
	@echo "  CC    test_aoa_triangulation"
	$(CC) $(CFLAGS) $< -o $@ -L$(LIB_DIR) -lmqp_full $(LDLIBS)

# Multi-module tests (noise, beamforming, CFAR, state machine, protocol)
test-modules: $(TEST_MODULES_BIN)
	@echo ""
	@echo "Running multi-module tests..."
	@echo "======================================"
	$(TEST_MODULES_BIN)

$(TEST_MODULES_BIN): $(TEST_MODULES_SRC) $(LIB_PATH) | $(BUILD_DIR)
	@echo "  CC    test_modules"
	$(CC) $(CFLAGS) $< -o $@ -L$(LIB_DIR) -lmqp_full $(LDLIBS)

# WPI campus integration tests
test-wpi: $(TEST_WPI_BIN)
	@echo ""
	@echo "Running WPI campus integration tests..."
	@echo "======================================"
	$(TEST_WPI_BIN)

$(TEST_WPI_BIN): $(TEST_WPI_SRC) $(LIB_PATH) | $(BUILD_DIR)
	@echo "  CC    test_wpi_integration"
	$(CC) $(CFLAGS) $< -o $@ -L$(LIB_DIR) -lmqp_full $(LDLIBS)

# Run all tests
test-all: test test-sigproc test-aoa test-modules test-wpi

# ============================================================================
# Module-specific Libraries
# ============================================================================

lib-music: $(LIB_DIR)/libmqp_music.a
$(LIB_DIR)/libmqp_music.a: $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(MUSIC_SOURCES)) | $(LIB_DIR)
	$(AR) rcs $@ $^

lib-signal: $(LIB_DIR)/libmqp_signal.a
$(LIB_DIR)/libmqp_signal.a: $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SIGNAL_SOURCES)) | $(LIB_DIR)
	$(AR) rcs $@ $^

lib-protocol: $(LIB_DIR)/libmqp_protocol.a
$(LIB_DIR)/libmqp_protocol.a: $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(PROTOCOL_SOURCES)) | $(LIB_DIR)
	$(AR) rcs $@ $^

lib-localization: $(LIB_DIR)/libmqp_localization.a
$(LIB_DIR)/libmqp_localization.a: $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(LOCALIZATION_SOURCES)) | $(LIB_DIR)
	$(AR) rcs $@ $^

# ============================================================================
# Installation
# ============================================================================

install: $(LIB_PATH)
	@echo "Installing library to $(INSTALL_LIB_DIR)..."
	install -d $(INSTALL_LIB_DIR)
	install -m 644 $(LIB_PATH) $(INSTALL_LIB_DIR)/
	@echo "Installing headers to $(INSTALL_INC_DIR)..."
	install -d $(INSTALL_INC_DIR)
	install -m 644 $(INC_DIR)/*.h $(INSTALL_INC_DIR)/
	@echo "Installation complete."

uninstall:
	rm -f $(INSTALL_LIB_DIR)/$(LIB_NAME)
	rm -rf $(INSTALL_INC_DIR)

# ============================================================================
# Maintenance
# ============================================================================

clean:
	@echo "Cleaning build artifacts..."
	rm -rf $(BUILD_DIR) $(LIB_DIR)

check-deps:
	@echo "Checking dependencies..."
	@command -v $(CC) >/dev/null 2>&1 || { echo "ERROR: GCC not found"; exit 1; }
	@pkg-config --exists openblas 2>/dev/null || echo "WARNING: OpenBLAS not found via pkg-config"
	@pkg-config --exists fftw3 2>/dev/null || echo "WARNING: FFTW3 not found via pkg-config"
	@pkg-config --exists libiio 2>/dev/null || echo "WARNING: LibIIO not found via pkg-config"
	@echo "Dependency check complete."

info:
	@echo "MQP Build Configuration"
	@echo "========================"
	@echo "Compiler: $(CC)"
	@$(CC) --version 2>/dev/null | head -1 || true
	@echo ""
	@echo "CFLAGS:  $(CFLAGS)"
	@echo "LDLIBS:  $(LDLIBS)"
	@echo ""
	@echo "Source files: $(words $(ALL_SOURCES))"
	@echo "  MUSIC/CFAR:         $(words $(MUSIC_SOURCES))"
	@echo "  Signal Processing:  $(words $(SIGNAL_SOURCES))"
	@echo "  Protocol Detection: $(words $(PROTOCOL_SOURCES))"
	@echo "  Localization:       $(words $(LOCALIZATION_SOURCES))"
	@echo "  SDR:                $(words $(SDR_SOURCES))"
	@echo ""
	@echo "Target: Raspberry Pi 5B (ARM Cortex-A76)"

help:
	@echo "MQP Full System Build"
	@echo "====================="
	@echo ""
	@echo "Targets:"
	@echo "  make all             Build library + app + tests"
	@echo "  make app             Build mqp_localize application"
	@echo "  make run             Build and run application (WPI debug mode)"
	@echo "  make lib             Build static library only"
	@echo "  make test            Build and run MUSIC unit tests"
	@echo "  make test-all        Run all tests"
	@echo "  make test-sdr        Build and run SDR tests (requires hardware)"
	@echo "  make clean           Remove build artifacts"
	@echo "  make install         Install library/headers (requires sudo)"
	@echo "  make check-deps      Verify dependencies"
	@echo "  make info            Display build configuration"
	@echo ""
	@echo "  make lib-music       Build MUSIC-only library"
	@echo "  make lib-signal      Build signal processing library"
	@echo "  make lib-protocol    Build protocol detection library"
	@echo "  make lib-localization Build localization library"
	@echo ""
	@echo "Options:"
	@echo "  DEBUG=1              Debug build (verbose, profiling, -O0)"
	@echo "  USE_SDR=1            Include ADALM-PLUTO SDR support"
	@echo ""
	@echo ""
	@echo "Application:"
	@echo "  ./build/mqp_localize -h               # Show app help"
	@echo "  ./build/mqp_localize -s wpi -b 0 -d   # WPI LTE Band 13 debug"
	@echo "  ./build/mqp_localize -s wpi -b 4       # WPI P25 Worcester PD"
	@echo ""
	@echo "Examples:"
	@echo "  make clean && make all               # Clean rebuild"
	@echo "  make DEBUG=1 test                    # Debug + test"
	@echo "  make USE_SDR=1 all                   # Build with SDR"
	@echo "  make USE_SDR=1 DEBUG=1 test-sdr      # Debug SDR test"
	@echo "  make USE_SDR=1 app                   # Build app with SDR"

# ============================================================================
# Advanced Targets
# ============================================================================

# Build with profiling
profile: CFLAGS += -pg -DMUSIC_PROFILE=1
profile: clean all

# Build with sanitizers
sanitize: CFLAGS += -fsanitize=address,undefined -fno-omit-frame-pointer
sanitize: LDLIBS += -fsanitize=address,undefined
sanitize: clean all

# ============================================================================
# Dependencies Installation (Raspberry Pi OS / Ubuntu ARM64):
#   sudo apt-get update
#   sudo apt-get install -y build-essential libopenblas-dev libfftw3-dev
#   sudo apt-get install -y libiio-dev libiio-utils   # For SDR support
#
# Performance targets:
#   MUSIC AOA:          <60ms
#   2D CFAR:            <50ms
#   Signal Processing:  <20ms
#   Protocol Detection: <30ms
#   Total Pipeline:     <200ms (real-time at 5 Hz)
# ============================================================================
