# mk/common.mk - Common definitions and platform detection

# Platform detection
UNAME_S := $(shell uname -s)

# Compiler and base flags
CC ?= cc
CFLAGS ?= -O3
CFLAGS += -Wall -Wextra

# Debug build support (append, don't replace)
ifeq ($(DEBUG), 1)
    CFLAGS := -O0 -g -DDEBUG -Wall -Wextra
endif

# Sanitizer support
ifeq ($(SANITIZE), 1)
    CFLAGS += -fsanitize=address,undefined
    LDFLAGS += -fsanitize=address,undefined
endif

# Directory structure
INCLUDE_DIR := include
SRC_DIR := src
EXAMPLES_DIR := examples
ASSETS_DIR := assets

# Library source and object
LIB_SRC := $(SRC_DIR)/b3d.c
LIB_OBJ := $(SRC_DIR)/b3d.o

# Library dependencies
LIB_DEPS := $(SRC_DIR)/math-toolkit.h $(INCLUDE_DIR)/b3d.h

# Common include path
INCLUDES := -I$(INCLUDE_DIR)

# Common libraries
LIBS := -lm

# Utility function: check if command exists
define has
$(shell command -v $(1) > /dev/null 2>&1 && echo 1)
endef

# Verbosity control (V=1 for verbose)
ifeq ($(V), 1)
    Q :=
    VECHO := @true
else
    Q := @
    VECHO := @echo
endif
