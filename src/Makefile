# Compiler selection: use CXX=g++ or CXX=clang++ (default: g++)
# Production build: AVX-512 IFMA. For integer fallback (debugging), use Makefile.fallback.
CXX ?= g++
INCLUDES = -I.. $(BENCHMARK_INC)

# Suppress warnings from BLST/vector code: static decls without defs, ignored attributes on __m512i, unknown pragmas (clang-only)
CXXFLAGS_WARN_SUPPRESS = -Wno-unused-function -Wno-ignored-attributes -Wno-unknown-pragmas
# Compiler-specific flags
ifneq ($(findstring clang++,$(CXX)),)
    CXXFLAGS = -std=c++20 -Wall -Wextra -Wno-sign-conversion $(CXXFLAGS_WARN_SUPPRESS) -mavx512ifma
    CXXFLAGS_BENCH = -std=c++20 -O3 -march=native -Wall -Wextra -Wno-sign-conversion $(CXXFLAGS_WARN_SUPPRESS) -mavx512ifma
else
    CXXFLAGS = -std=c++2a -Wall -Wextra -Wno-sign-conversion $(CXXFLAGS_WARN_SUPPRESS) -mavx512ifma
    CXXFLAGS_BENCH = -std=c++2a -O3 -march=native -Wall -Wextra -Wno-sign-conversion $(CXXFLAGS_WARN_SUPPRESS) -mavx512ifma
endif

BLST_LIBS = -lgmp -lgmpxx ../blst/libblst_small.a
BENCHMARK_INC ?=
BENCHMARK_LIBS ?= -lbenchmark -lpthread

TARGET_REDUCE = test_elementwise_reduce
SOURCE_REDUCE = test_elementwise_reduce.cpp
TARGET_INVERSION = test_inversion
SOURCE_INVERSION = test_inversion.cpp
TARGET_TEST_PAIRING = test_pairing
SOURCE_TEST_PAIRING = test_pairing.cpp
TARGET_TEST_EC = test_ec
SOURCE_TEST_EC = test_ec.cpp
TARGET_TEST_BLS12_381_FIELD_MUL = test_bls12_381_field_mul
SOURCE_TEST_BLS12_381_FIELD_MUL = test_bls12_381_field_mul.cpp
TARGET_TEST_BLS12_381_PAIRING = test_bls12_381_pairing
SOURCE_TEST_BLS12_381_PAIRING = test_bls12_381_pairing.cpp
TARGET_TEST_EC_BN254 = test_ec_bn254
SOURCE_TEST_EC_BN254 = test_ec_bn254.cpp
TARGET_TEST_SCALAR_MULT = test_scalar_mult
SOURCE_TEST_SCALAR_MULT = test_scalar_mult.cpp
TARGET_BENCH_EC_BN254 = bench_ec_bn254
SOURCE_BENCH_EC_BN254 = bench_ec_bn254.cpp
TARGET_BENCH_BLS12_381 = bench_bls12_381
SOURCE_BENCH_BLS12_381 = bench_bls12_381.cpp bench_bls12_381_impl.hpp bench_bls12_381_common.hpp

HEADERS = bounds.hpp constexpr_utils.hpp bounded_element.hpp ring_element.hpp \
          elementwise_reduce.hpp change_base.hpp bounded_ring.hpp \
          fp12.hpp final_exponentiation.hpp miller.hpp fp2.hpp inversion.hpp \
          ec.hpp conversion_inversion.hpp scalar_mult.hpp

.PHONY: all clean test test-reduce test-inversion test-pairing test-ec \
        test-bls12-381-field-mul test-bls12-381-pairing test-ec-bn254 test-scalar-mult \
        bench-ec-bn254 bench-bls12-381 run

all: $(TARGET_REDUCE) $(TARGET_INVERSION) $(TARGET_TEST_EC) $(TARGET_TEST_SCALAR_MULT)

$(TARGET_REDUCE): $(SOURCE_REDUCE) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SOURCE_REDUCE) -o $(TARGET_REDUCE)

$(TARGET_INVERSION): $(SOURCE_INVERSION) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SOURCE_INVERSION) -o $(TARGET_INVERSION) -lgmp -lgmpxx

$(TARGET_TEST_PAIRING): $(SOURCE_TEST_PAIRING) $(HEADERS)
	$(CXX) $(CXXFLAGS) -D__ADX__ $(INCLUDES) $(SOURCE_TEST_PAIRING) -o $(TARGET_TEST_PAIRING) $(BLST_LIBS)

$(TARGET_TEST_EC): $(SOURCE_TEST_EC) $(HEADERS)
	$(CXX) $(CXXFLAGS) -D__ADX__ $(INCLUDES) $(SOURCE_TEST_EC) -o $(TARGET_TEST_EC) $(BLST_LIBS)

$(TARGET_TEST_BLS12_381_FIELD_MUL): $(SOURCE_TEST_BLS12_381_FIELD_MUL) $(HEADERS)
	$(CXX) $(CXXFLAGS) -D__ADX__ $(INCLUDES) $(SOURCE_TEST_BLS12_381_FIELD_MUL) -o $(TARGET_TEST_BLS12_381_FIELD_MUL) -lgmp -lgmpxx

$(TARGET_TEST_BLS12_381_PAIRING): $(SOURCE_TEST_BLS12_381_PAIRING) pairing.hpp bls12_381_ring_params.hpp
	$(CXX) $(CXXFLAGS) -D__ADX__ $(INCLUDES) $(SOURCE_TEST_BLS12_381_PAIRING) -o $(TARGET_TEST_BLS12_381_PAIRING) $(BLST_LIBS)

$(TARGET_TEST_EC_BN254): $(SOURCE_TEST_EC_BN254) $(HEADERS)
	$(CXX) $(CXXFLAGS) -D__ADX__ $(INCLUDES) $(SOURCE_TEST_EC_BN254) -o $(TARGET_TEST_EC_BN254) $(BLST_LIBS)

$(TARGET_TEST_SCALAR_MULT): $(SOURCE_TEST_SCALAR_MULT) $(HEADERS)
	$(CXX) $(CXXFLAGS) -D__ADX__ $(INCLUDES) $(SOURCE_TEST_SCALAR_MULT) -o $(TARGET_TEST_SCALAR_MULT) $(BLST_LIBS)

$(TARGET_BENCH_BLS12_381): $(SOURCE_BENCH_BLS12_381) $(HEADERS) bls12_381_ring_params.hpp
	$(CXX) $(CXXFLAGS_BENCH) -D__ADX__ $(INCLUDES) bench_bls12_381.cpp -o $(TARGET_BENCH_BLS12_381) $(BLST_LIBS) $(BENCHMARK_LIBS)

$(TARGET_BENCH_EC_BN254): $(SOURCE_BENCH_EC_BN254) $(HEADERS)
	$(CXX) $(CXXFLAGS_BENCH) -D__ADX__ $(INCLUDES) $(SOURCE_BENCH_EC_BN254) -o $(TARGET_BENCH_EC_BN254) -lgmp -lgmpxx $(BENCHMARK_LIBS)

test: $(TARGET_REDUCE) $(TARGET_INVERSION)
	./$(TARGET_REDUCE)
	./$(TARGET_INVERSION)

test-reduce: $(TARGET_REDUCE)
	./$(TARGET_REDUCE)

test-inversion: $(TARGET_INVERSION)
	./$(TARGET_INVERSION)

test-pairing: $(TARGET_TEST_PAIRING)
	./$(TARGET_TEST_PAIRING)

test-ec: $(TARGET_TEST_EC)
	./$(TARGET_TEST_EC)

test-bls12-381-field-mul: $(TARGET_TEST_BLS12_381_FIELD_MUL)
	./$(TARGET_TEST_BLS12_381_FIELD_MUL)

test-bls12-381-pairing: $(TARGET_TEST_BLS12_381_PAIRING)
	./$(TARGET_TEST_BLS12_381_PAIRING)

test-ec-bn254: $(TARGET_TEST_EC_BN254)
	./$(TARGET_TEST_EC_BN254)

test-scalar-mult: $(TARGET_TEST_SCALAR_MULT)
	./$(TARGET_TEST_SCALAR_MULT)

bench-ec-bn254: $(TARGET_BENCH_EC_BN254)
	./$(TARGET_BENCH_EC_BN254)

bench-bls12-381: $(TARGET_BENCH_BLS12_381)
	./$(TARGET_BENCH_BLS12_381)

run: test

clean:
	rm -f $(TARGET_REDUCE) $(TARGET_INVERSION) \
	      $(TARGET_TEST_PAIRING) $(TARGET_TEST_EC) \
	      $(TARGET_TEST_BLS12_381_FIELD_MUL) $(TARGET_TEST_BLS12_381_PAIRING) \
	      $(TARGET_TEST_EC_BN254) $(TARGET_TEST_SCALAR_MULT) \
	      $(TARGET_BENCH_EC_BN254) $(TARGET_BENCH_BLS12_381)
