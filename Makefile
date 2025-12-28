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
	$(Q)rm -f $(ALL_EXAMPLES_CLEAN) $(LIB_OBJ) $(VOXEL_OBJ) $(TESTS) $(BENCHMARKS) tests/test-gen $(MATH_GEN_TEST_H)

# Clean everything including generated source files
distclean: cleanall
	$(Q)rm -f src/math-gen.inc examples/globe-data.h

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

# Check examples use b3d-math.h wrappers (not raw math.h functions)
check-math-usage:
	$(VECHO) "  CHECK\tVerifying examples use b3d-math.h"
	$(Q)python3 scripts/check-math-usage.py

# Run unit tests
check: check-math-usage $(TESTS)
	$(VECHO) "  TEST\tRunning unit tests"
	$(Q)set -e; for t in $(TESTS); do $$t; done
	@echo ""
	@echo "All tests passed"

# Run benchmarks
bench: $(BENCHMARKS)
	$(VECHO) "  BENCH\tRunning performance benchmarks"
	$(Q)set -e; for b in $(BENCHMARKS); do $$b; done

# Run all tests (unit + benchmarks)
test-all: check bench

# Code generation for dual float/fixed-point math
generate: $(MATH_GEN_H)

$(MATH_GEN_H): $(MATH_DSL) $(MATH_GEN)
	$(VECHO) "  GEN\t$@"
	$(Q)python3 $(MATH_GEN) --dsl $(MATH_DSL) -o $@

# Test version with _gen suffix to avoid conflicts with existing functions
$(MATH_GEN_TEST_H): $(MATH_DSL) $(MATH_GEN)
	$(VECHO) "  GEN\t$@ (test)"
	$(Q)python3 $(MATH_GEN) --dsl $(MATH_DSL) --suffix _gen -o $@

# Verify generated code compiles and produces correct results
tests/test-gen: tests/test-gen.c $(MATH_GEN_TEST_H) $(LIB_DEPS)
	$(VECHO) "  CC\t$@"
	$(Q)$(CC) $(CFLAGS) -DB3D_FLOAT_POINT $(INCLUDES) -Isrc $< -o $@ $(LIBS)

check-gen: tests/test-gen
	$(VECHO) "  TEST\tRunning generated code tests"
	$(Q)./tests/test-gen

# Snapshot generation (outputs to examples/*.png)
update-snapshots: $(ALL_EXAMPLES)
	$(VECHO) "  GEN\tGenerating snapshots"
	$(Q)./scripts/gen-snapshots.sh $(SDL2_EXAMPLES_ALL)

# Download Natural Earth data and create binary asset
# Requires: pip3 install pyshp requests
update-globe-data:
	$(VECHO) "  GET\t$(GLOBE_BINARY)"
	$(Q)python3 scripts/gen-globe-data.py download -o $(GLOBE_BINARY)

.PHONY: all clean cleanall rebuild config check check-math-usage bench test-all generate check-gen update-snapshots update-globe-data $(BUILD_TARGETS) $(RUN_TARGETS)
