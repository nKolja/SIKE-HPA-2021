/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: benchmarking/testing isogeny-based key encapsulation mechanism
*********************************************************************************************/ 


// Benchmark and test parameters  
#if defined(OPTIMIZED_GENERIC_IMPLEMENTATION) || (TARGET == TARGET_ARM) 
    #define BENCH_LOOPS        5      // Number of iterations per bench 
    #define TEST_LOOPS         5      // Number of iterations per test
#else
    #define BENCH_LOOPS       100       
    #define TEST_LOOPS        10      
#endif

int main(int argc, char* argv[])
{
    int Status = PASSED;
    Status = hamming_weight_algorithm(argc, argv);             // Key extraction algorithm
    if (Status != PASSED) {
        printf("\n\n   Error detected! \n\n");
        return FAILED;
    }


    return Status;
}