# Minimal working examples for VROOM RNS Montgomery multiplication.
#
# Default: AVX-512 IFMA (same as src/Makefile). For machines without IFMA, use
#   make FALLBACK=1
#
# Requires: g++/clang++ with C++20, GMP (libgmp, libgmpxx).

CXX ?= g++
INCLUDES = -I..
CXXFLAGS_WARN_SUPPRESS = -Wno-unused-function -Wno-ignored-attributes -Wno-unknown-pragmas

ifeq ($(FALLBACK),1)
  VECTOR_FLAGS = -DUSE_AVX512_FALLBACK
else
  VECTOR_FLAGS = -mavx512ifma
endif

ifneq ($(findstring clang++,$(CXX)),)
  CXXFLAGS = -std=c++20 -O2 -Wall -Wextra -Wno-sign-conversion $(CXXFLAGS_WARN_SUPPRESS) $(VECTOR_FLAGS)
else
  CXXFLAGS = -std=c++2a -O2 -Wall -Wextra -Wno-sign-conversion $(CXXFLAGS_WARN_SUPPRESS) $(VECTOR_FLAGS)
endif

LIBS = -lgmp -lgmpxx

EXAMPLES = 01_parameter_setup 02_singular_modmul 03_sum_of_products

.PHONY: all clean $(EXAMPLES)

all: $(EXAMPLES)

01_parameter_setup: 01_parameter_setup.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< -o $@ $(LIBS)

02_singular_modmul: 02_singular_modmul.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< -o $@ $(LIBS)

03_sum_of_products: 03_sum_of_products.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< -o $@ $(LIBS)

clean:
	rm -f $(EXAMPLES)
