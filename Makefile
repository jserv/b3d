# Makefile for bootleg3D
# Builds the library and all example programs

CC = cc
CFLAGS = -O3 -Wall -Wextra

# SDL2 configuration (auto-detect via sdl2-config)
SDL_CFLAGS = $(shell sdl2-config --cflags)
SDL_LIBS = $(shell sdl2-config --libs)

# Directory structure
INCLUDE_DIR = include
SRC_DIR = src
EXAMPLES_DIR = examples

# Library source and object
LIB_SRC = $(SRC_DIR)/b3d.c
LIB_OBJ = $(SRC_DIR)/b3d.o

# Example source files
CUBES_SRC = $(EXAMPLES_DIR)/cubes.c
OBJ_SRC = $(EXAMPLES_DIR)/obj.c
FPS_SRC = $(EXAMPLES_DIR)/fps.c
TERRAIN_SRC = $(EXAMPLES_DIR)/terrain.c
ASCII_SRC = $(EXAMPLES_DIR)/ascii.c
DONUT_SRC = $(EXAMPLES_DIR)/donut.c
LENA3D_SRC = $(EXAMPLES_DIR)/lena3d.c

# Executable names
CUBES_EXE = cubes
OBJ_EXE = obj
FPS_EXE = fps
TERRAIN_EXE = terrain
ASCII_EXE = ascii
DONUT_EXE = donut
LENA3D_EXE = lena3d

# All targets
ALL_TARGETS = $(CUBES_EXE) $(OBJ_EXE) $(FPS_EXE) $(TERRAIN_EXE) $(ASCII_EXE) $(DONUT_EXE) $(LENA3D_EXE)

# Default target - build all examples
all: $(ALL_TARGETS)

# Build library object file
$(LIB_OBJ): $(LIB_SRC) $(SRC_DIR)/math.h $(INCLUDE_DIR)/b3d.h
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) -c $< -o $@

# Examples that require SDL2
$(CUBES_EXE): $(CUBES_SRC) $(LIB_OBJ)
	$(CC) $(CFLAGS) $(SDL_CFLAGS) -I$(INCLUDE_DIR) $^ -o $@ $(SDL_LIBS) -lm

$(OBJ_EXE): $(OBJ_SRC) $(LIB_OBJ)
	$(CC) $(CFLAGS) $(SDL_CFLAGS) -I$(INCLUDE_DIR) $^ -o $@ $(SDL_LIBS) -lm

$(FPS_EXE): $(FPS_SRC) $(LIB_OBJ)
	$(CC) $(CFLAGS) $(SDL_CFLAGS) -I$(INCLUDE_DIR) $^ -o $@ $(SDL_LIBS) -lm

$(TERRAIN_EXE): $(TERRAIN_SRC) $(LIB_OBJ)
	$(CC) $(CFLAGS) $(SDL_CFLAGS) -I$(INCLUDE_DIR) $^ -o $@ $(SDL_LIBS) -lm

$(ASCII_EXE): $(ASCII_SRC) $(LIB_OBJ)
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) $^ -o $@ -lm

$(DONUT_EXE): $(DONUT_SRC) $(LIB_OBJ)
	$(CC) $(CFLAGS) $(SDL_CFLAGS) -I$(INCLUDE_DIR) $^ -o $@ $(SDL_LIBS) -lm

$(LENA3D_EXE): $(LENA3D_SRC) $(LIB_OBJ)
	$(CC) $(CFLAGS) $(SDL_CFLAGS) -I$(INCLUDE_DIR) $^ -o $@ $(SDL_LIBS) -lm

# Aliases for building individual examples
build-cubes: $(CUBES_EXE)
build-obj: $(OBJ_EXE)
build-fps: $(FPS_EXE)
build-terrain: $(TERRAIN_EXE)
build-ascii: $(ASCII_EXE)
build-donut: $(DONUT_EXE)
build-lena3d: $(LENA3D_EXE)

# Clean build artifacts
clean:
	rm -f $(ALL_TARGETS) $(LIB_OBJ)

# Clean everything including generated assets
cleanall: clean
	rm -f assets/cube.png

# Rebuild everything
rebuild: clean all

# Run examples from top-level directory
run-cubes: $(CUBES_EXE)
	./$(CUBES_EXE)

run-obj: $(OBJ_EXE)
	./$(OBJ_EXE)

run-fps: $(FPS_EXE)
	./$(FPS_EXE)

run-terrain: $(TERRAIN_EXE)
	./$(TERRAIN_EXE)

run-ascii: $(ASCII_EXE)
	./$(ASCII_EXE)

run-donut: $(DONUT_EXE)
	./$(DONUT_EXE)

run-lena3d: $(LENA3D_EXE)
	./$(LENA3D_EXE)

# List all targets
list:
	@echo "Available targets:"
	@echo "  all         - Build all examples"
	@echo "  build-cubes - Build cubes example (requires SDL2)"
	@echo "  build-obj   - Build obj example (requires SDL2)"
	@echo "  build-fps   - Build fps example (requires SDL2)"
	@echo "  build-terrain - Build heightmap example (requires SDL2)"
	@echo "  build-ascii - Build ASCII cube example (no SDL2)"
	@echo "  build-donut - Build torus demo (requires SDL2)"
	@echo "  build-lena3d - Build Lena heightfield demo (requires SDL2)"
	@echo "  run-cubes  - Run cubes example"
	@echo "  run-obj    - Run obj example"
	@echo "  run-fps    - Run fps example"
	@echo "  run-terrain - Run heightmap example"
	@echo "  run-ascii - Run ASCII cube example (prints to stdout)"
	@echo "  run-donut - Run torus demo"
	@echo "  run-lena3d - Run Lena heightfield demo"
	@echo "  clean      - Remove all built executables and object files"
	@echo "  cleanall   - Remove executables, object files, and generated assets"
	@echo "  rebuild    - Clean and rebuild all examples"

.PHONY: all clean cleanall rebuild list build-cubes build-obj build-fps build-terrain build-ascii build-donut build-lena3d run-cubes run-obj run-fps run-terrain run-ascii run-donut run-lena3d
