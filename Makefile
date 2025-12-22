# Configuration options:
#   ENABLE_SDL2=0/1 - Toggle SDL2 examples (auto-detected)
#   DEBUG=1         - Debug build with symbols
#   SANITIZE=1      - Enable address/undefined sanitizers
#   V=1             - Verbose output

# Default target (must be before includes that define targets)
.DEFAULT_GOAL := all

# Tests
TESTS := tests/math-fixed tests/math-float

# Include modular build components
include mk/common.mk
include mk/sdl2.mk
include mk/examples.mk

# Default target
all: $(ALL_EXAMPLES)

# Library object file
$(LIB_OBJ): $(LIB_SRC) $(LIB_DEPS)
	$(VECHO) "  CC\t$@"
	$(Q)$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Clean build artifacts (all possible examples, regardless of config)
clean:
	$(VECHO) "  CLEAN"
	$(Q)rm -f $(ALL_EXAMPLES_CLEAN) $(LIB_OBJ) $(TESTS)

# Clean everything including generated assets
cleanall: clean
	$(VECHO) "  CLEANALL"
	$(Q)rm -f $(ASSETS_DIR)/*.png

# Rebuild everything
rebuild: clean all

# Show configuration
config:
	@echo "B3D build configuration:"
	@echo "  Platform:    $(UNAME_S)"
	@echo "  Compiler:    $(CC)"
	@echo "  CFLAGS:      $(CFLAGS)"
	@echo "  LDFLAGS:     $(LDFLAGS)"
	@echo "  SDL2:        $(if $(filter 1,$(ENABLE_SDL2)),enabled,disabled)"
	@echo "  Examples:    $(ALL_EXAMPLES)"

# Tests
tests/math-fixed: tests/test-math.c $(LIB_DEPS)
	$(VECHO) "  CC\t$@"
	$(Q)$(CC) $(CFLAGS) $(INCLUDES) -Isrc $< -o $@ $(LIBS)

tests/math-float: tests/test-math.c $(LIB_DEPS)
	$(VECHO) "  CC\t$@ (float)"
	$(Q)$(CC) $(CFLAGS) -DB3D_FLOAT_POINT $(INCLUDES) -Isrc $< -o $@ $(LIBS)

check: $(TESTS)
	$(VECHO) "  RUN\t$(TESTS)"
	$(Q)set -e; for t in $(TESTS); do $$t; done

.PHONY: all clean cleanall rebuild config test $(BUILD_TARGETS) $(RUN_TARGETS)
