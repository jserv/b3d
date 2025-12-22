# Configuration options:
#   ENABLE_SDL2=0/1 - Toggle SDL2 examples (auto-detected)
#   DEBUG=1         - Debug build with symbols
#   SANITIZE=1      - Enable address/undefined sanitizers
#   V=1             - Verbose output

# Default target (must be before includes that define targets)
.DEFAULT_GOAL := all

# Tests (unit tests)
TESTS := tests/math-fixed tests/math-float tests/test-api
# Benchmarks (performance tests)
BENCHMARKS := tests/test-perf

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
	$(Q)rm -f $(ALL_EXAMPLES_CLEAN) $(LIB_OBJ) $(TESTS) $(BENCHMARKS)

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

tests/test-api: tests/test-api.c $(LIB_DEPS) $(LIB_OBJ)
	$(VECHO) "  CC\t$@"
	$(Q)$(CC) $(CFLAGS) $(INCLUDES) $< $(LIB_OBJ) -o $@ $(LIBS)

tests/test-perf: tests/test-perf.c $(LIB_DEPS) $(LIB_OBJ)
	$(VECHO) "  CC\t$@"
	$(Q)$(CC) $(CFLAGS) $(INCLUDES) $< $(LIB_OBJ) -o $@ $(LIBS)

# Run unit tests
check: $(TESTS)
	$(VECHO) "  TEST\tRunning unit tests"
	$(Q)set -e; for t in $(TESTS); do $$t; done
	@echo ""
	@echo "All unit tests passed"

# Run benchmarks
bench: $(BENCHMARKS)
	$(VECHO) "  BENCH\tRunning performance benchmarks"
	$(Q)set -e; for b in $(BENCHMARKS); do $$b; done

# Run all tests (unit + benchmarks)
test-all: check bench

.PHONY: all clean cleanall rebuild config check bench test-all $(BUILD_TARGETS) $(RUN_TARGETS)
