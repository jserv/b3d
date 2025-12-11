# mk/sdl2.mk - SDL2 detection and configuration

# Only run detection if ENABLE_SDL2 not explicitly set
ifeq ($(origin ENABLE_SDL2), undefined)
    SDL2_CONFIG := $(call has, sdl2-config)
    ifeq ($(SDL2_CONFIG), 1)
        ENABLE_SDL2 := 1
    else
        ENABLE_SDL2 := 0
        $(warning SDL2 not found via sdl2-config, SDL2 examples disabled)
    endif
endif

# Configure SDL2 flags when enabled
ifeq ($(ENABLE_SDL2), 1)
    SDL_CFLAGS := $(shell sdl2-config --cflags 2>/dev/null)
    SDL_LIBS := $(shell sdl2-config --libs 2>/dev/null)
else
    SDL_CFLAGS :=
    SDL_LIBS :=
endif
