# Convenience wrapper around CMake — NOT the build system itself.
# The real build system is CMake (see CMakeLists.txt); these targets just
# invoke it with the common flag combinations so you don't have to remember
# them. Everything below delegates to `cmake`; there are no hand-written
# compile rules.
#
#   make              optimized (Release) build
#   make help         list all targets
#   make test         build + run the test suite
#   make install      install to PREFIX (default: CMake's default prefix)
#
# Overridable variables:
#   BUILD=build       build directory
#   PREFIX=/usr/local install prefix (for `make install`)
#   JOBS=8            parallel build jobs (default: detected core count)

BUILD  ?= build
JOBS   ?= $(shell getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)
CMAKE  ?= cmake

# Common invocations
CONFIGURE = $(CMAKE) -B $(BUILD)
BUILDCMD  = $(CMAKE) --build $(BUILD) -j $(JOBS)

.DEFAULT_GOAL := release

.PHONY: help release debug dev test benchmark fuzz shared native install clean

help:              ## list available targets
	@echo "carquet — CMake convenience wrapper"
	@echo
	@grep -E '^[a-z][a-zA-Z]*:.*##' $(MAKEFILE_LIST) \
		| sed -E 's/^([a-zA-Z]+):.*## (.*)/\1|\2/' \
		| awk -F'|' '{printf "  \033[1m%-12s\033[0m %s\n", $$1, $$2}'
	@echo
	@echo "Variables: BUILD=$(BUILD) PREFIX=<prefix> JOBS=$(JOBS)"

release:           ## optimized build (default)
	$(CONFIGURE) -DCMAKE_BUILD_TYPE=Release
	$(BUILDCMD)

debug:             ## debug build (assertions, no optimization)
	$(CONFIGURE) -DCMAKE_BUILD_TYPE=Debug
	$(BUILDCMD)

dev:               ## build all dev targets (tests, examples, benchmarks, interop)
	$(CONFIGURE) -DCMAKE_BUILD_TYPE=Debug -DCARQUET_BUILD_DEV=ON
	$(BUILDCMD)

test:              ## build tests and run ctest
	$(CONFIGURE) -DCMAKE_BUILD_TYPE=Debug -DCARQUET_BUILD_TESTS=ON
	$(BUILDCMD)
	ctest --test-dir $(BUILD) --output-on-failure

benchmark:         ## build benchmark programs
	$(CONFIGURE) -DCMAKE_BUILD_TYPE=Release -DCARQUET_BUILD_BENCHMARKS=ON
	$(BUILDCMD)

fuzz:              ## build fuzz targets
	$(CONFIGURE) -DCMAKE_BUILD_TYPE=Debug -DCARQUET_BUILD_FUZZ=ON
	$(BUILDCMD)

shared:            ## build shared library
	$(CONFIGURE) -DCMAKE_BUILD_TYPE=Release -DCARQUET_BUILD_SHARED=ON
	$(BUILDCMD)

native:            ## optimized build tuned for this host (-march=native)
	$(CONFIGURE) -DCMAKE_BUILD_TYPE=Release -DCARQUET_NATIVE_ARCH=ON
	$(BUILDCMD)

install: release   ## install to PREFIX (e.g. make install PREFIX=/usr/local)
	$(CMAKE) --install $(BUILD)$(if $(PREFIX), --prefix $(PREFIX),)

clean:             ## remove the build directory
	rm -rf $(BUILD)
