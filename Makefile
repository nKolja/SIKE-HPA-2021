####  Makefile for compilation on Unix-like operative systems. Updated for with side-channel analysis functions.  ####

OPT=-O3     # Optimization option by default

CC=clang
ifeq "$(CC)" "gcc"
    COMPILER=gcc
else ifeq "$(CC)" "clang"
    COMPILER=clang
endif

ARCHITECTURE=_AMD64_
USE_OPT_LEVEL=_FAST_
ifeq "$(ARCH)" "x64"
    ARCHITECTURE=_AMD64_
    USE_OPT_LEVEL=_FAST_
else ifeq "$(ARCH)" "x86"
    ARCHITECTURE=_X86_
    USE_OPT_LEVEL=_GENERIC_
else ifeq "$(ARCH)" "s390x"
    ARCHITECTURE=_S390X_
    USE_OPT_LEVEL=_GENERIC_
else ifeq "$(ARCH)" "ARM"
    ARCHITECTURE=_ARM_
    USE_OPT_LEVEL=_GENERIC_
    ARM_TARGET=YES
else ifeq "$(ARCH)" "ARM64"
    ARCHITECTURE=_ARM64_
    USE_OPT_LEVEL=_FAST_
    ARM_TARGET=YES
endif

ifeq "$(OPT_LEVEL)" "GENERIC"
    USE_OPT_LEVEL=_GENERIC_
endif

ifeq "$(ARM_TARGET)" "YES"
    ARM_SETTING=-lrt
endif

ifeq "$(ARCHITECTURE)" "_AMD64_"
    ifeq "$(USE_OPT_LEVEL)" "_FAST_"
        MULX=-D _MULX_
        ifeq "$(USE_MULX)" "FALSE"
            MULX=
        else
            ADX=-D _ADX_
            ifeq "$(USE_ADX)" "FALSE"
                ADX=
            endif
        endif
    endif
endif

AR=ar rcs
RANLIB=ranlib

ADDITIONAL_SETTINGS=-march=native
ifeq "$(CC)" "clang"
ifeq "$(ARM_TARGET)" "YES"
    ADDITIONAL_SETTINGS=
endif
endif
ifeq "$(ARCHITECTURE)" "_S390X_"
	ADDITIONAL_SETTINGS=-march=z10
endif

CFLAGS=$(OPT) -std=gnu11 -I/usr/local/include $(ADDITIONAL_SETTINGS) -D $(ARCHITECTURE) -D __NIX__ -D $(USE_OPT_LEVEL) $(MULX) $(ADX)
LDFLAGS=-lm -lgmp
ifeq "$(USE_OPT_LEVEL)" "_GENERIC_"
    EXTRA_OBJECTS_434=objs434/fp_generic.o
else ifeq "$(USE_OPT_LEVEL)" "_FAST_"
ifeq "$(ARCHITECTURE)" "_AMD64_"
    EXTRA_OBJECTS_434=objs434/fp_x64.o objs434/fp_x64_asm.o
else ifeq "$(ARCHITECTURE)" "_ARM64_"
    EXTRA_OBJECTS_434=objs434/fp_arm64.o objs434/fp_arm64_asm.o
endif
endif
OBJECTS_434=objs434/P434.o $(EXTRA_OBJECTS_434) objs/random.o objs/fips202.o

all: lib434 tests KATS Bob_Key Alice_Key make_points

objs434/%.o: src/P434/%.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@

ifeq "$(USE_OPT_LEVEL)" "_GENERIC_"		
    objs434/fp_generic.o: src/P434/generic/fp_generic.c
	    $(CC) -c $(CFLAGS) src/P434/generic/fp_generic.c -o objs434/fp_generic.o
else ifeq "$(USE_OPT_LEVEL)" "_FAST_"
ifeq "$(ARCHITECTURE)" "_AMD64_"		
    objs434/fp_x64.o: src/P434/AMD64/fp_x64.c
	    $(CC) -c $(CFLAGS) src/P434/AMD64/fp_x64.c -o objs434/fp_x64.o

    objs434/fp_x64_asm.o: src/P434/AMD64/fp_x64_asm.S
	    $(CC) -c $(CFLAGS) src/P434/AMD64/fp_x64_asm.S -o objs434/fp_x64_asm.o
else ifeq "$(ARCHITECTURE)" "_ARM64_"	
    objs434/fp_arm64.o: src/P434/ARM64/fp_arm64.c
	    $(CC) -c $(CFLAGS) src/P434/ARM64/fp_arm64.c -o objs434/fp_arm64.o

    objs434/fp_arm64_asm.o: src/P434/ARM64/fp_arm64_asm.S
	    $(CC) -c $(CFLAGS) src/P434/ARM64/fp_arm64_asm.S -o objs434/fp_arm64_asm.o
endif
endif

objs/random.o: src/random/random.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) src/random/random.c -o objs/random.o

objs/fips202.o: src/sha3/fips202.c
	$(CC) -c $(CFLAGS) src/sha3/fips202.c -o objs/fips202.o

lib434: $(OBJECTS_434)
	rm -rf lib434 sike434 sidh434
	mkdir lib434 sike434 sidh434
	$(AR) lib434/libsidh.a $^
	$(RANLIB) lib434/libsidh.a

tests: lib434 #USED TO COMPUTE HAMMING WEIGHTS OF ELLIPTIC CURVE POINTS USED IN SIDE-CHANNEL ATTACK
	$(CC) $(CFLAGS) -L./lib434 tests/test_p434_hamming_weights.c tests/test_extras.c -lsidh $(LDFLAGS) -o hamming_weight_computation $(ARM_SETTING)

Bob_Key: lib434 #USED TO CREATE A KEYPAIR FOR BOB
	$(CC) $(CFLAGS) -L./lib434 tests/test_p434_bob.c tests/test_extras.c -lsidh $(LDFLAGS) -o bob_key $(ARM_SETTING)

Alice_Key: lib434 #USED TO CREATE A KEYPAIR FOR ALICE
	$(CC) $(CFLAGS) -L./lib434 tests/test_p434_alice.c tests/test_extras.c -lsidh $(LDFLAGS) -o alice_key $(ARM_SETTING)

make_points: lib434 #USED TO CREATE POINTS FOR THE MONTGOMERY LADDER AND TO UPDATE THEM AFTER GUESSING BITS
	$(CC) $(CFLAGS) -L./lib434 tests/test_p434_make_points.c tests/test_extras.c -lsidh $(LDFLAGS) -o make_starting_points $(ARM_SETTING)
# AES
AES_OBJS=objs/aes.o objs/aes_c.o

objs/%.o: tests/aes/%.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@

lib434_for_KATs: $(OBJECTS_434) $(AES_OBJS)
	$(AR) lib434/libsidh_for_testing.a $^
	$(RANLIB) lib434/libsidh_for_testing.a

KATS: lib434_for_KATs 
	$(CC) $(CFLAGS) -L./lib434 tests/PQCtestKAT_kem434.c tests/rng/rng.c -lsidh_for_testing $(LDFLAGS) -o sike434/PQCtestKAT_kem $(ARM_SETTING)

check: tests

.PHONY: clean

clean:
	rm -rf *.req objs434* objs503* objs610* objs751* objs lib434* lib503* lib610* lib751* sidh434* sidh503* sidh610* sidh751* sike434* sike503* sike610* sike751* arith_tests-* hamming_weight_computation weights/* points/* bob_key alice_key make_starting_points

