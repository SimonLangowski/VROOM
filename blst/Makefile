# Makefile for building blst subset (multi_scalar.c and pairing.c)
# Only x86-64 assembly is needed
# Uses separate compilation units for better optimization

# Allow CC to be overridden from command line, otherwise auto-detect
# Supports both C (gcc/clang) and C++ (g++/clang++) compilers
ifeq ($(origin CC),default)
  CC := $(shell command -v clang 2>/dev/null || command -v gcc 2>/dev/null || echo gcc)
endif

# Detect if C++ compiler (ends with ++)
IS_CPP := $(findstring ++,$(CC))

AR = ar
RANLIB = ranlib

# Base flags - match original build.sh flags exactly
CFLAGS = -O2 -fno-builtin -fPIC -Wall -Wextra -I. -D__x86_64__

# Check for ADX support
HAS_ADX := $(shell grep -q "adx" /proc/cpuinfo 2>/dev/null && echo yes || echo no)
ifeq ($(HAS_ADX),yes)
  CFLAGS += -D__ADX__
endif

# Compiler-specific flags
ifneq ($(IS_CPP),)
  # C++ compiler flags - compile as C++
  ifeq ($(findstring clang++,$(CC)),clang++)
    CFLAGS += -std=c++17 -Wno-implicit-function-declaration
    AR = llvm-ar
    RANLIB = llvm-ranlib
  else
    # g++ flags - compile as C++
    CFLAGS += -std=c++17
  endif
else
  # C compiler flags
  ifeq ($(findstring clang,$(CC)),clang)
    # Clang-specific flags: allow implicit function declarations (like gcc does)
    CFLAGS += -Wno-implicit-function-declaration
    AR = llvm-ar
    RANLIB = llvm-ranlib
  else
    # GCC-specific flags
    CFLAGS += -std=gnu99
  endif
endif

# Source files - compile as single unit like original build.sh
SERVER_SRC = server_small.c
SERVER_OBJ = server_small.o

# Assembly file - single unit like original build/assembly.S
ASSEMBLY_SRC = assembly_small.S
ASSEMBLY_OBJ = assembly_small.o

# Benchmark program (optional - only if source exists)
BENCHMARK_SRC = benchmark.c
ifneq ($(wildcard $(BENCHMARK_SRC)),)
  BENCHMARK_OBJ = benchmark.o
  BENCHMARK = benchmark
  BENCHMARK_TARGETS = $(BENCHMARK)
else
  BENCHMARK_OBJ =
  BENCHMARK =
  BENCHMARK_TARGETS =
endif

# Empirical basis program (C++)
EMPIRICAL_SRC = empirical_basis.cpp
EMPIRICAL_OBJ = empirical_basis.o
EMPIRICAL = empirical_basis

# Test program for final exponentiation
TEST_FINAL_EXP_SRC = test_final_exponentiation.cpp
TEST_FINAL_EXP_OBJ = test_final_exponentiation.o
TEST_FINAL_EXP = test_final_exponentiation

# Test program for inversion parts
TEST_INVERSION_PARTS_SRC = test_inversion_parts.cpp
TEST_INVERSION_PARTS_OBJ = test_inversion_parts.o
TEST_INVERSION_PARTS = test_inversion_parts

# FP12 test program (optional - only if source exists)
FP12_SRC = fp12.cpp
ifneq ($(wildcard $(FP12_SRC)),)
  FP12_OBJ = fp12.o
  FP12 = fp12
  FP12_TARGETS = $(FP12)
else
  FP12_OBJ =
  FP12 =
  FP12_TARGETS =
endif

# FP12Ring test generator
TEST_FP12RING_GENERATE_SRC = test_fp12ring_generate.cpp
TEST_FP12RING_GENERATE_OBJ = test_fp12ring_generate.o
TEST_FP12RING_GENERATE = test_fp12ring_generate

# Benchmark pairing program (optional - only if source exists)
BENCH_PAIRING_BLST_SRC = bench_pairing_blst.cpp
ifneq ($(wildcard $(BENCH_PAIRING_BLST_SRC)),)
  BENCH_PAIRING_BLST_OBJ = bench_pairing_blst.o
  BENCH_PAIRING_BLST = bench_pairing_blst
  BENCH_PAIRING_BLST_TARGETS = $(BENCH_PAIRING_BLST)
else
  BENCH_PAIRING_BLST_OBJ =
  BENCH_PAIRING_BLST =
  BENCH_PAIRING_BLST_TARGETS =
endif

# Google Benchmark harness over blst API (bindings/blst.h)
BENCH_BLST_COMPLETE_SRC = bench_blst_complete.cpp
BENCH_BLST_GLUE_SRC = bench_blst_glue.c
ifneq ($(wildcard $(BENCH_BLST_COMPLETE_SRC)),)
  BENCH_BLST_COMPLETE_OBJ = bench_blst_complete.o
  BENCH_BLST_GLUE_OBJ = bench_blst_glue.o
  BENCH_BLST_COMPLETE = bench_blst_complete
  BENCH_BLST_COMPLETE_TARGETS = $(BENCH_BLST_COMPLETE)
else
  BENCH_BLST_COMPLETE_OBJ =
  BENCH_BLST_GLUE_OBJ =
  BENCH_BLST_COMPLETE =
  BENCH_BLST_COMPLETE_TARGETS =
endif

# Output library
LIBRARY = libblst_small.a

# All targets: library + only programs whose sources exist
all: $(LIBRARY) $(EMPIRICAL) $(TEST_FINAL_EXP) $(TEST_FP12RING_GENERATE) $(TEST_INVERSION_PARTS) $(BENCHMARK_TARGETS) $(FP12_TARGETS) $(BENCH_PAIRING_BLST_TARGETS) $(BENCH_BLST_COMPLETE_TARGETS)

$(LIBRARY): $(SERVER_OBJ) $(ASSEMBLY_OBJ)
	$(AR) rc $(LIBRARY) $(SERVER_OBJ) $(ASSEMBLY_OBJ)
	$(RANLIB) $(LIBRARY)

# Compile server_small.c as single unit (includes all .c files)
# For C++ compilers, compile as C++ - no -x c flag
# For C compilers, compile as C
$(SERVER_OBJ): $(SERVER_SRC)
ifneq ($(IS_CPP),)
	$(CC) $(CFLAGS) -x c++ -c $(SERVER_SRC) -o $(SERVER_OBJ)
else
	$(CC) $(CFLAGS) -c $(SERVER_SRC) -o $(SERVER_OBJ)
endif

# Compile assembly_small.S as single unit (includes all .s files)
# Assembly must be compiled as assembly, but needs C preprocessor for #include
# Remove language-specific flags for assembly
ASM_CFLAGS = $(filter-out -x c -x c++ -std=gnu99 -Wno-implicit-function-declaration,$(CFLAGS))
$(ASSEMBLY_OBJ): $(ASSEMBLY_SRC)
	$(CC) -x assembler-with-cpp $(ASM_CFLAGS) -c $(ASSEMBLY_SRC) -o $(ASSEMBLY_OBJ)

# Compile benchmark.c - use same language mode as server (only if source exists)
ifneq ($(BENCHMARK_TARGETS),)
$(BENCHMARK_OBJ): $(BENCHMARK_SRC)
ifneq ($(IS_CPP),)
	$(CC) $(CFLAGS) -x c++ -c $(BENCHMARK_SRC) -o $(BENCHMARK_OBJ)
else
	$(CC) $(CFLAGS) -c $(BENCHMARK_SRC) -o $(BENCHMARK_OBJ)
endif

$(BENCHMARK): $(BENCHMARK_OBJ) $(LIBRARY)
	$(CC) $(CFLAGS) -o $(BENCHMARK) $(BENCHMARK_OBJ) -L. -lblst_small
endif

# Compile empirical_basis.cpp - always compile as C++
# Filter out C-specific flags and add C++ flags
# Suppress unused-function warnings to see actual errors
CXXFLAGS = $(filter-out -std=gnu99,$(CFLAGS)) -std=c++17 -Wno-unused-function
$(EMPIRICAL_OBJ): $(EMPIRICAL_SRC)
	g++ $(CXXFLAGS) -c $(EMPIRICAL_SRC) -o $(EMPIRICAL_OBJ)

$(EMPIRICAL): $(EMPIRICAL_OBJ) $(LIBRARY)
	g++ $(CXXFLAGS) -o $(EMPIRICAL) $(EMPIRICAL_OBJ) -L. -lblst_small

$(TEST_FINAL_EXP_OBJ): $(TEST_FINAL_EXP_SRC)
	g++ $(CXXFLAGS) -c $(TEST_FINAL_EXP_SRC) -o $(TEST_FINAL_EXP_OBJ)

$(TEST_FINAL_EXP): $(TEST_FINAL_EXP_OBJ) $(LIBRARY)
	g++ $(CXXFLAGS) -o $(TEST_FINAL_EXP) $(TEST_FINAL_EXP_OBJ) -L. -lblst_small

$(TEST_INVERSION_PARTS_OBJ): $(TEST_INVERSION_PARTS_SRC)
	g++ $(CXXFLAGS) -c $(TEST_INVERSION_PARTS_SRC) -o $(TEST_INVERSION_PARTS_OBJ)

$(TEST_INVERSION_PARTS): $(TEST_INVERSION_PARTS_OBJ) $(LIBRARY)
	g++ $(CXXFLAGS) -o $(TEST_INVERSION_PARTS) $(TEST_INVERSION_PARTS_OBJ) -L. -lblst_small

# FP12 and bench_pairing_blst (only if sources exist)
ifneq ($(FP12_TARGETS),)
$(FP12_OBJ): $(FP12_SRC)
	g++ $(CXXFLAGS) -c $(FP12_SRC) -o $(FP12_OBJ)

$(FP12): $(FP12_OBJ) $(LIBRARY)
	g++ $(CXXFLAGS) -o $(FP12) $(FP12_OBJ) -L. -lblst_small
endif

$(TEST_FP12RING_GENERATE_OBJ): $(TEST_FP12RING_GENERATE_SRC)
	g++ $(CXXFLAGS) -c $(TEST_FP12RING_GENERATE_SRC) -o $(TEST_FP12RING_GENERATE_OBJ)

$(TEST_FP12RING_GENERATE): $(TEST_FP12RING_GENERATE_OBJ) $(LIBRARY)
	g++ $(CXXFLAGS) -o $(TEST_FP12RING_GENERATE) $(TEST_FP12RING_GENERATE_OBJ) -L. -lblst_small

ifneq ($(BENCH_PAIRING_BLST_TARGETS),)
$(BENCH_PAIRING_BLST_OBJ): $(BENCH_PAIRING_BLST_SRC)
	g++ $(CXXFLAGS) -c $(BENCH_PAIRING_BLST_SRC) -o $(BENCH_PAIRING_BLST_OBJ)

$(BENCH_PAIRING_BLST): $(BENCH_PAIRING_BLST_OBJ) $(LIBRARY)
	g++ $(CXXFLAGS) -o $(BENCH_PAIRING_BLST) $(BENCH_PAIRING_BLST_OBJ) -L. -lblst_small -lbenchmark -lpthread -lgmp -lgmpxx
endif

ifneq ($(BENCH_BLST_COMPLETE_TARGETS),)
.PHONY: gbench
# C glue: libblst_small inlines Fp ops; bench calls blst_fp_* from bindings/blst.h
GLUE_CFLAGS = -O2 -fno-builtin -fPIC -Wall -Wextra -I. -D__x86_64__ -std=gnu99
ifeq ($(HAS_ADX),yes)
  GLUE_CFLAGS += -D__ADX__
endif

$(BENCH_BLST_GLUE_OBJ): $(BENCH_BLST_GLUE_SRC) fields.h vect.h consts.h
	gcc $(GLUE_CFLAGS) -c $(BENCH_BLST_GLUE_SRC) -o $(BENCH_BLST_GLUE_OBJ)

$(BENCH_BLST_COMPLETE_OBJ): $(BENCH_BLST_COMPLETE_SRC) bindings/blst.h
	g++ $(CXXFLAGS) -c $(BENCH_BLST_COMPLETE_SRC) -o $(BENCH_BLST_COMPLETE_OBJ)

$(BENCH_BLST_COMPLETE): $(BENCH_BLST_COMPLETE_OBJ) $(BENCH_BLST_GLUE_OBJ) $(LIBRARY)
	g++ $(CXXFLAGS) -o $(BENCH_BLST_COMPLETE) $(BENCH_BLST_COMPLETE_OBJ) $(BENCH_BLST_GLUE_OBJ) -L. -lblst_small -lbenchmark -lpthread

gbench: $(BENCH_BLST_COMPLETE)
endif

clean:
	rm -f $(SERVER_OBJ) $(ASSEMBLY_OBJ) $(LIBRARY) $(BENCHMARK_OBJ) $(BENCHMARK) $(EMPIRICAL_OBJ) $(EMPIRICAL) $(TEST_FINAL_EXP_OBJ) $(TEST_FINAL_EXP) $(FP12_OBJ) $(FP12) $(TEST_FP12RING_GENERATE_OBJ) $(TEST_FP12RING_GENERATE) $(BENCH_PAIRING_BLST_OBJ) $(BENCH_PAIRING_BLST) $(BENCH_BLST_COMPLETE_OBJ) $(BENCH_BLST_GLUE_OBJ) $(BENCH_BLST_COMPLETE) $(TEST_INVERSION_PARTS_OBJ) $(TEST_INVERSION_PARTS)

.PHONY: all clean
