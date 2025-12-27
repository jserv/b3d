# mk/examples.mk - Example build rules

# Explicitly listed examples (add new examples here)
SDL2_EXAMPLES_ALL := ball cubes donut fps gears lena3d lighting obj terrain voxel

# Enable SDL2 examples only if SDL2 available
ifeq ($(ENABLE_SDL2), 1)
    SDL2_EXAMPLES := $(SDL2_EXAMPLES_ALL)
else
    SDL2_EXAMPLES :=
endif

# All examples to build
ALL_EXAMPLES := $(SDL2_EXAMPLES)

# All possible examples (for clean target)
ALL_EXAMPLES_CLEAN := $(SDL2_EXAMPLES_ALL)

# Auto-generate source mappings: <name>_SRC
$(foreach ex,$(ALL_EXAMPLES_CLEAN),$(eval $(ex)_SRC := $(EXAMPLES_DIR)/$(ex).c))

# Extra dependencies (header-only libs, assets, etc.)
obj_DEPS := $(INCLUDE_DIR)/b3d-obj.h
fps_DEPS := $(INCLUDE_DIR)/b3d-obj.h
voxel_DEPS := $(INCLUDE_DIR)/b3d-voxel.h $(INCLUDE_DIR)/b3d-obj.h

# Voxel module object file
VOXEL_OBJ := $(SRC_DIR)/b3d-voxel.o
VOXEL_SRC := $(SRC_DIR)/b3d-voxel.c

# Build rule for voxel object
$(VOXEL_OBJ): $(VOXEL_SRC) $(INCLUDE_DIR)/b3d-voxel.h $(INCLUDE_DIR)/b3d.h
	$(VECHO) "  CC\t$@"
	$(Q)$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Per-example extra objects (beyond LIB_OBJ)
voxel_OBJS := $(VOXEL_OBJ)

# Per-example CFLAGS (if needed)
# fps_CFLAGS := -DFPS_DEBUG

# Per-example LDFLAGS (if needed)
# obj_LDFLAGS := -lpng

# Pattern rule for SDL2 examples
define sdl2_example_rule
$(1): $$($(1)_SRC) $$(LIB_OBJ) $$($(1)_OBJS) $$($(1)_DEPS)
	$$(VECHO) "  CC\t$$@"
	$$(Q)$$(CC) $$(CFLAGS) $$($(1)_CFLAGS) $$(SDL_CFLAGS) $$(INCLUDES) $$< $$(LIB_OBJ) $$($(1)_OBJS) -o $$@ $$(SDL_LIBS) $$(LIBS) $$(LDFLAGS) $$($(1)_LDFLAGS)
endef

# Generate build rules for all examples
$(foreach ex,$(SDL2_EXAMPLES),$(eval $(call sdl2_example_rule,$(ex))))

# Build aliases: build-<example>
$(foreach ex,$(ALL_EXAMPLES),$(eval build-$(ex): $(ex)))

# Run targets: run-<example>
define run_rule
run-$(1): $(1)
	$$(VECHO) "  RUN\t$(1)"
	$$(Q)./$(1)
endef
$(foreach ex,$(ALL_EXAMPLES),$(eval $(call run_rule,$(ex))))

# Phony declarations for build/run targets
BUILD_TARGETS := $(addprefix build-,$(ALL_EXAMPLES))
RUN_TARGETS := $(addprefix run-,$(ALL_EXAMPLES))
