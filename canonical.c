#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "fasta.h"

// =================================================================================================
//     Helper Functions
// =================================================================================================

// print a k-mer in binary format
void printkmer(u_int64_t u)
{
    for (int i = 0; i < 64; ++i) {
        if (i > 0 && i % 2 == 0) {
            printf(" ");
        }
        u_int64_t const pos = 1ULL << (64 - i - 1);
        if (u & pos) {
            printf("1");
        } else {
            printf("0");
        }
    }
    printf("\n");
}

/**
 * @brief Calculate the power `base^exp` for positive integer values.
 *
 * This is considerably faster than pow(), and does not have any rounding or converion issues.
 * However, it overflows quite easily. The function does not check whether the desired power
 * actually fits within `u_int64_t`.
 */
u_int64_t int_pow(u_int64_t base, u_int8_t exp)
{
    // Using Exponentiation by squaring, see
    // http://stackoverflow.com/a/101613/4184258
    u_int64_t result = 1;
    while (exp) {
        if (exp & 1) {
            result *= base;
        }
        exp >>= 1;
        base *= base;
    }
    return result;
}

/**
 * @brief Compute the total number of possible k-mers for a given @p k.
 */
u_int64_t number_of_kmers(u_int8_t k)
{
    if (k == 0 || k >= 32) {
        fprintf(stderr, "ERROR: Can only compute number of k-mers for k in [1,32).");
        return 0;
    }
    u_int64_t n = 1;
    for (u_int64_t i = 0; i < k; ++i) {
        n *= 4;
    }
    return n;
}

/**
 * @brief Compute the number of canonical k-mers for a given k and nucleotide alphabet.
 */
u_int64_t number_of_canonical_kmers(u_int8_t k)
{
    // We need distinct approaches for even and odd values, due to palindromes.
    // It might be easier to just have a hard coded table... but this way is more approachable.
    if (k == 0 || k > 32) {
        fprintf(stderr, "ERROR: Can only compute number of canonical k-mers for k in [1,32].");
    } else if (k % 2 == 0) {
        // Even numbers, need to add palindromes.
        // We use base 2 here, and instead of dividing the result by 2 in the end, we subtract 1
        // from the exponent, in order to avoid overflowing for the case k=32.
        // The original (overflowing) equation from the paper is commented out below for reference.
        return int_pow(2, 2 * k - 1) + int_pow(2, 2 * k / 2 - 1);
        // return ( int_pow( 4, k ) + int_pow( 4, k / 2 )) / 2;
    } else {
        // Odd numbers. No overflow for the valid range.
        return int_pow(4, k) / 2;
    }
    return 0;
}

/**
 * @brief Compute the number of palindromes (under reverse complmenet) that exist
 * for a given @p k and nucleotide alphabet.
 *
 * This is `0` for odd values of @p k, and `4^(k/2)` for even values of @p k.
 */
u_int64_t number_of_palindromes(u_int8_t k)
{
    // Edge and special cases.
    if (k == 0 || k > 32) {
        fprintf(stderr, "ERROR: Can only compute number of palindromes for k in [1,32].");
        return 0;
    } else if (k % 2 != 0) {
        // No palindromes for odd k
        return 0;
    }
    return int_pow(4, k / 2);
}

/**
 * @brief Compute the reverse complement of a give @p kmer.
 */
u_int64_t reverse_complement(u_int64_t kmer, u_int8_t k)
{
    // Adapted from Kraken2 at https://github.com/DerrickWood/kraken2/blob/master/src/mmscanner.cc
    // which itself adapted this for 64-bit DNA use from public domain code at
    // https://graphics.stanford.edu/~seander/bithacks.html#ReverseParallel

    // Reverse bits (leaving bit pairs intact, as those represent nucleotides):
    // Swap consecutive pairs, then nibbles, then bytes, then byte pairs, then halves of 64-bit word
    u_int64_t value = kmer;
    value = ((value & 0xCCCCCCCCCCCCCCCCUL) >> 2) | ((value & 0x3333333333333333UL) << 2);
    value = ((value & 0xF0F0F0F0F0F0F0F0UL) >> 4) | ((value & 0x0F0F0F0F0F0F0F0FUL) << 4);
    value = ((value & 0xFF00FF00FF00FF00UL) >> 8) | ((value & 0x00FF00FF00FF00FFUL) << 8);
    value = ((value & 0xFFFF0000FFFF0000UL) >> 16) | ((value & 0x0000FFFF0000FFFFUL) << 16);
    value = (value >> 32) | (value << 32);

    // Finally, complement, and shift to correct position, removing the invalid lower bits.
    int const bitwidth = sizeof(u_int64_t) * 8;
    value = ((~value) >> (bitwidth - 2 * k));
    return value;
}

// =================================================================================================
//     Minimal Canonical Encoding Functions
// =================================================================================================

typedef struct {
    int k;

    u_int64_t allones;
    int bitwidth;

    // Precomputed powers of 4
    u_int64_t four_to_the_k_half_plus_one;
    u_int64_t twice_four_to_the_k_half;
    u_int64_t gap_shift;

    // max + 1, since max = 8 * sizeof(u_int64_t)
    u_int64_t zeromasks[65];
    u_int64_t remaindermasks[34]; // max k + 2
} Bitmasks;

// Replace markers for the function R for each type of specifying pair.
// We code those as a lookup table, where each entry is a single word
// containing the four bits of the following list in their LSBs.
// We store those as the type of our underlying data, so that we can
// directly shift those values to the position where they are needed.
// *  0 A..A            -> 0110
// *  1 A..C            -> 0101
// *  2 A..G            -> 0100
// #  3 A..T palindrome -> 0000
// *  4 C..A            -> 1000
// *  5 C..C            -> 0111
// #  6 C..G palindrome -> 0000
// #  7 C..T -> A..G    -> 0100
// *  8 G..A            -> 1001
// #  9 G..C palindrome -> 0000
// # 10 G..G -> C..C    -> 0111
// # 11 G..T -> A..C    -> 0101
// # 12 T..A palindrome -> 0000
// # 13 T..C -> G..A    -> 1001
// # 14 T..G -> C..A    -> 1000
// # 15 T..T -> A..A    -> 0110
static const u_int64_t replace[16] = {
    0x06, 0x05, 0x04, 0x00, 0x08, 0x07, 0x00, 0x04, 0x09, 0x00, 0x07, 0x05, 0x00, 0x09, 0x08, 0x06
};

// Markers to check if we need to encode the forward or the reverse complement.
static const char reverse[16] = { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1 };

void initialize_bitmasks(Bitmasks* bm, int k)
{
    bm->k = k;
    bm->bitwidth = 8 * sizeof(u_int64_t);

    // Precompute constants for correcting the gap sizes in the index
    bm->four_to_the_k_half_plus_one = int_pow(4, k / 2 + 1);
    bm->twice_four_to_the_k_half = 2 * int_pow(4, k / 2);
    bm->gap_shift = (2 * ((k + 1) / 2) - 1);

    // Precompute masks
    bm->allones = (((1ULL << 32) - 1) << 32) + ((1ULL << 32) - 1);

    // Initialize zeromasks
    // 0 .. 0 1 .. 1
    bm->zeromasks[0] = bm->allones;
    for (int i = 1; i <= bm->bitwidth; i++) {
        bm->zeromasks[i] = bm->zeromasks[i - 1] >> 1;
    }

    // After we have identified the specifying pair of characters, we need to extract
    // the remainder, see encode_prime(). We here precompute a mask to do that.
    // For instance, for k==7, the relevant entries are shaped like this:
    //
    //     remainder_mask_[2] == 00 .. 00 11 11 11 11 00
    //     remainder_mask_[4] == 00 .. 00 00 11 11 00 00
    //     remainder_mask_[6] == 00 .. 00 00 11 00 00 00
    //
    // We only ever need to access entries at even indices, as this is indexed per bit,
    // and we use index access to the starting bit of the characters.
    // Lastly, as explained in encode_prime(), we also might access entries
    // beyond the given triangle of 1s, so we fill those with zeros here.
    bm->remaindermasks[0] = bm->allones;
    for (size_t i = 1; i <= k; i++) {
        u_int64_t const zeromask = bm->allones >> (bm->bitwidth - 2 * k + i);
        u_int64_t const onemask = bm->allones << i;
        bm->remaindermasks[i] = zeromask & onemask;
    }
    for (size_t i = k + 1; i < 34; ++i) {
        bm->remaindermasks[i] = 0;
    }
}

// Compute encoding where only setting the bits according to specifying case and
// subtracting gaps is missing, i.e., enc prime.
u_int64_t encode_prime(Bitmasks const* bm, u_int64_t kmer, int l)
{
    int k = bm->k;

    // This uses a mask of the form 0..01..1 (l trailing ones), to extract
    // the relevant bits on the right, and invert (complement) them.
    // u_int64_t const zeromask = (l == 0 ? 0 : bm->allones >> (bm->bitwidth - l));
    u_int64_t const zeromask = bm->zeromasks[bm->bitwidth - l];
    u_int64_t const right = (kmer & zeromask) ^ zeromask;

    // No remainder left? We could just return here, but in our tests, the introduced
    // branching is more expensive than unconditionally executing the below bit operations,
    // so we have deactivated this check here. Recommended to be tested on your hardware.
    // if( l + 2 >= bm->k ) {
    //     return right;
    // }

    // Assert that the values are as expected.
    // assert(l <= bm->k);
    // assert(l % 2 == 0);

    // Use the remainder mask (consisting of ones in the middle) to extract the bits
    // in between the specifying pair, then shift the remainder to the correct position.
    // The mask contains 0 after index k, so that if we have l+2 >= k (no remainder),
    // we just get a zero here, which does nothing to our result.
    u_int64_t const remainder = (kmer & bm->remaindermasks[l + 2]) >> 2;
    return right | remainder;
}

// compute enc_r_c
u_int64_t encode(Bitmasks const* bm, u_int64_t kmer, u_int64_t rckmer)
{
    // hash
    u_int64_t kmercode;
    int k = bm->k;

    // Get the length of the symmetric prefix/suffix, in num of characters, i.e., 2x num of bits.
    // Then, l is the bit index of the char that is the specifying case for the k-mer.
    // Calling ctz(0) is undefined, which is the case for palindromes; we catch this
    // below by checking for sym==0, as this is faster than an additional check here.
    u_int64_t sym = kmer ^ rckmer;
    int l = __builtin_ctzll(sym) / 2 * 2;

    if (sym == 0) {
        // Palindrome -> nothing to do. Can only occurr in even k.
        // assert(k % 2 == 0);

        // We use l = k here, as in a palindrome, l will overshoot due to the ctl call,
        // so we limit it here to the range that we are interested in.
        l = k;
        kmercode = encode_prime(bm, kmer, l);
    } else if (l < k - 1) {
        // Not just single character in the middle, i.e., we have a specifying pair.

        // There are 16 possible combinations of two characters from ACGT.
        // We here extract the first two asymmetric characters (the specifying pair, i.e. 2x2 bits)
        // to build a pattern for a lookup of which combination we have in the kmer.
        // This is done by shifting the relevant bits of the pair to the LSBs of the pattern.
        unsigned char pattern = 0;
        pattern |= (kmer >> (2 * k - l - 4)) & 0x0C;
        pattern |= (kmer >> l) & 0x03;
        // assert(pattern < 16);

        // Check which case we need for the initial hash, based on the pattern we found.
        if (reverse[pattern]) {
            kmercode = encode_prime(bm, rckmer, l);
        } else {
            kmercode = encode_prime(bm, kmer, l);
        }

        // Set positions l+1, l+2, l+3 and l+4 according to the specifying pair pattern,
        // which is called R in the manuscript.
        // Similar to above, we can avoid any branching here by directly shifting
        // the replace mask bits to the needed positions. If the replace mask is 0 for the
        // given pattern, we shift a zero, which just does nothing.
        kmercode |= (replace[pattern] << (2 * k - l - 4));
    } else if (l == k - 1) {
        // Single character in the middle. Can only occurr in odd k.
        // assert(k % 2 == 1);

        // We are interested in the bits at the central character,
        // which (given that we have l == k - 1 here) are located at:
        //     2*k - l - 1 == k
        //     2*k - l - 2 == k - 1
        // Use these bits to encode A/T -> 0 and C/G -> 1:
        //     A = 00 -> 0
        //     C = 01 -> 1
        //     G = 10 -> 1
        //     T = 11 -> 0
        // Depending on the combination of those two bits, we want to set a bit in kmercode.
        // In particular, we want to set the same bit as the second of the two above positions,
        // but only if both bit positions are different (C or G). For this, we first obtain
        // both bits of the kmer, and use XOR to see if they are different. To this end,
        // bit1 is shifted by 1 so that it is in the same positon as bit2.
        // The result of this XOR is a single bit indicating if we have C/G or A/T at the position,
        // and it is already in the correct position to be set in kmercode.
        kmercode = encode_prime(bm, kmer, l);
        u_int64_t const bit1 = (kmer & (1ULL << (k))) >> 1;
        u_int64_t const bit2 = (kmer & (1ULL << (k - 1)));
        kmercode |= bit1 ^ bit2;
    } else {
        // This branch is taken only if l >= k, which however implies that the kmer
        // is a palindrome, which is already checked as the first branch.
        assert(false);
    }

    // subtract gaps
    // 2*(k//2-l-1) ones followed by k-2 zeros
    if (l <= k - 4) {
        u_int64_t gaps = bm->zeromasks[bm->bitwidth - (k / 2 * 2 - l - 2)];
        kmercode -= (gaps << bm->gap_shift);
    }

    // subtract gap in code due to specifying middle position
    if (k % 2 == 1 && kmercode >= bm->four_to_the_k_half_plus_one) {
        kmercode -= bm->twice_four_to_the_k_half;
    }

    return kmercode;
}

// =================================================================================================
//     Process: Minimal Canonical Encoding
// =================================================================================================

// extract all k-mers, encode by enc^r_c, assign to b bins
int process_string(char* s, int k, int* bins, int b)
{

    char warned = 0; //warn first time, an unsupported character is skipped
    u_int64_t kmer = 0;
    u_int64_t rckmer = 0;
    u_int64_t const pos1 = 1ull << (2 * k - 2);
    u_int64_t const pos2 = 1ull << (2 * k - 1);

    u_int64_t maxrank = int_pow(2, 2 * k - 1);
    if (k % 2 == 0) {
        // Even k -> need to consider palindromes.
        // This is not the actual max rank. It is chosen such that the bucket sizes are a bit smaller;
        // such that palindromic k-mers (which do not appear in pairs) contribute as much as the
        // other canonical k-mers
        // This is the same as above, but we state it explicitly here so indicate that
        // this can be adjusted if needed.
        // Using a power of 2 here to express (4^k) / 2 without overflow for k==32.
        maxrank = int_pow(2, 2 * k - 1);
    }

    Bitmasks bm;
    initialize_bitmasks(&bm, k);

    // go through the string, process each k-mer by shift/update, enc_r_c, and assign to bins
    for (int i = 0; i < strlen(s); i++) {

        //update kmer by current character
        kmer = kmer << 2;
        switch (s[i]) {
        case 'A':
            break;
        case 'C':
            kmer += 1;
            break;
        case 'G':
            kmer += 2;
            break;
        case 'T':
            kmer += 3;
            break;
        default:
            if (!warned) {
                fprintf(stderr, "Warning: unsupported character replaced by A: %c\n", s[i]);
                warned = 1;
            }
            break;
        }

        // clear unused bits
        kmer &= bm.zeromasks[bm.bitwidth - 2 * k];

        // update reverse complement kmer
        rckmer = rckmer >> 2;
        switch (s[i]) {
        case 'T':
            break;
        case 'G':
            rckmer |= pos1;
            break;
        case 'C':
            rckmer |= pos2;
            break;
        case 'A':
            rckmer |= pos1 | pos2;
            break;
        default:
            rckmer |= pos1 | pos2;
            break;
        }

        // Not yet done filling the first kmer
        if (i < k - 1) {
            continue;
        }

        // hash
        u_int64_t kmercode = encode(&bm, kmer, rckmer);

        // assign to bin
        u_int64_t p = b * kmercode / maxrank;

        if (p == b) {
            // palindromes comes not in pairs like the other canonical k-mers. In order to ensure a equal distribution, bucket sizes are chosen as if there were only half as many palindromes possible. This results in ranks that would correspond to a max+1st bin. These k-mers are assigned to bin 0 where the palindromes are "missing".
            bins[0]++;
        } else {
            bins[p]++;
        }
    }

    return 0;
}

// =================================================================================================
//     Process: Standard 2-bit Encoding
// =================================================================================================

// extract all k-mers, encode by standard 2-bit encoding, assign to b bins
int process_string_std(char* s, int k, int* bins, int b)
{

    u_int64_t kmer = 0;
    u_int64_t rckmer = 0;
    u_int64_t const rem = (k == 32 ? 0ull : 3ull) << (2 * k);
    u_int64_t const pos1 = 1ull << (2 * k - 2);
    u_int64_t const pos2 = 1ull << (2 * k - 1);
    u_int64_t const maxrank = int_pow(2, 2 * k);
    char warned = 0; //warn first time, an unsupported character is skipped

    for (int i = 0; i < strlen(s); i++) {

        //update kmer by current character
        kmer = kmer << 2;
        switch (s[i]) {
        case 'A':
            break;
        case 'C':
            kmer += 1;
            break;
        case 'G':
            kmer += 2;
            break;
        case 'T':
            kmer += 3;
            break;
        default:
            if (!warned) {
                fprintf(stderr, "Warning: unsupported character replaced by A: %c\n", s[i]);
                warned = 1;
            }
            break;
        }

        // clear unused bits
        kmer &= ~rem;

        // update reverse complement kmer
        rckmer = rckmer >> 2;
        switch (s[i]) {
        case 'T':
            break;
        case 'G':
            rckmer |= pos1;
            break;
        case 'C':
            rckmer |= pos2;
            break;
        case 'A':
            rckmer |= pos1 | pos2;
            break;
        default:
            rckmer |= pos1 | pos2;
            break;
        }

        if (i >= k - 1) {

            u_int64_t can = kmer < rckmer ? kmer : rckmer;

            // assign to bin
            u_int64_t p = b * can / maxrank;
            // int p = can % t ;
            bins[p]++;
        }
    }

    return 0;
}

// =================================================================================================
//     Tests
// =================================================================================================

// Function to generate a random 64-bit integer
u_int64_t random_64bit_int()
{
    u_int64_t random_number = 0;

    // Combine several calls to rand() to construct a 64-bit integer
    // rand() is usually a 15-bit random number generator or larger,
    // so we stack enough of those to get our bits.
    random_number ^= ((u_int64_t)rand() << 48);
    random_number ^= ((u_int64_t)rand() << 32);
    random_number ^= ((u_int64_t)rand() << 16);
    random_number ^= ((u_int64_t)rand());

    return random_number;
}

// Simple helper to print an error if the given value is not true.
// We could use normal assert(), but that could be deactivated if not compiled as debug.
void test_assert(bool value)
{
    if (!value) {
        perror("Failed test");
        exit(EXIT_FAILURE);
    }
}

// Function to test all kmers up to a given size
void test_all_small_kmers()
{
    Bitmasks bm;

    // Test several different lengths of kmers
    for (u_int64_t k = 1; k < 12; ++k) {
        u_int64_t num_canon_kmers = number_of_canonical_kmers(k);
        u_int64_t num_palindromes = number_of_palindromes(k);
        initialize_bitmasks(&bm, k);

        // Create an array to count canocical indices, initialized to 0
        u_int64_t* counts = (u_int64_t*)calloc(num_canon_kmers, sizeof(u_int64_t));
        if (counts == NULL) {
            perror("Failed to allocate memory for counts array");
            exit(EXIT_FAILURE);
        }

        // Test all kmers of that length
        for (u_int64_t i = 0; i < number_of_kmers(k); ++i) {
            // Create the kmer, by simply using i as the binary coding.
            // This way, each possilbe kmer will be encountered exactly once.
            u_int64_t kmer = i;
            u_int64_t rc = reverse_complement(kmer, k);

            // Get its index from the encoding.
            u_int64_t index = encode(&bm, kmer, rc);

            // The index needs to match the one of the reverse complement
            test_assert(index == encode(&bm, rc, kmer));

            // Increment the count of that index, checking that we are in bounds.
            test_assert(index < num_canon_kmers);
            ++counts[index];
        }

        // Test that all bins got the number of kmers that we expect.
        test_assert(num_canon_kmers == num_canon_kmers); // Counts size is num_canon_kmers
        u_int64_t cnt = 0;
        for (u_int64_t i = 0; i < num_canon_kmers; ++i) {
            // For palindromes: the first 4^(k/2)/2 entries are only set once.
            if (k % 2 == 0 && i < num_palindromes) {
                test_assert(counts[i] == 1);
            } else {
                test_assert(counts[i] == 2);
            }
            cnt += counts[i];
        }
        test_assert(cnt == number_of_kmers(k));

        // Free the counts array after use
        free(counts);
    }
}

void test_large_kmers()
{
    // Seed the random number generator
    srand(time(NULL));
    int k = 32;
    Bitmasks bm;
    initialize_bitmasks(&bm, k);

    // Test large sizes of k for the boundaries.
    // Here, we cannot enumerate all values, so we just test the basic canonical property.
    for (int i = 0; i < 100000; ++i) {
        // Make a random kmer
        u_int64_t kmer = random_64bit_int();
        u_int64_t rc = reverse_complement(kmer, k);

        // The index needs to match the one of the reverse complement
        test_assert(encode(&bm, kmer, rc) == encode(&bm, rc, kmer));
    }
}

void test_speed(int k)
{
    // Seed the random number generator
    srand(time(NULL));
    Bitmasks bm;
    initialize_bitmasks(&bm, k);

    // Define the number of random integers we want to generate
    const int NUM_RANDOM_INTS = 100000000;

    // Allocate memory for arrays of kmers and their reverse complements.
    u_int64_t* kmers = (u_int64_t*)malloc(NUM_RANDOM_INTS * sizeof(u_int64_t));
    u_int64_t* rcs = (u_int64_t*)malloc(NUM_RANDOM_INTS * sizeof(u_int64_t));
    if (kmers == NULL || rcs == NULL) {
        perror("Failed to allocate memory for kmers or rcs");
        exit(EXIT_FAILURE);
    }

    // Generate random 64-bit integers and store them in the array
    for (int i = 0; i < NUM_RANDOM_INTS; ++i) {
        kmers[i] = random_64bit_int() >> (64 - 2 * k);
        rcs[i] = reverse_complement(kmers[i], k);
    }

    // Start high-resolution timer
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // Test that the encoding is the same for the kmer and its rc.
    // That's our speed test, hence encoding twice the number of kmers of the array.
    for (int i = 0; i < NUM_RANDOM_INTS; ++i) {
        test_assert(encode(&bm, kmers[i], rcs[i]) == encode(&bm, rcs[i], kmers[i]));
    }

    // Calculate the elapsed time in seconds, and the number of encodings per sec we achieved.
    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    unsigned long enc_per_sec = 2 * NUM_RANDOM_INTS / elapsed_time;
    printf("k==%i, time: %fs, enc/s: %lu\n", k, elapsed_time, enc_per_sec);

    // Free the allocated memory
    free(kmers);
    free(rcs);
}

// =================================================================================================
//     Main
// =================================================================================================

int main(int argc, char* argv[])
{
    // Run the test cases instead of the main processing.
    // test_all_small_kmers();
    // test_large_kmers();
    // test_speed(15);
    // test_speed(16);
    // return 0;

    // default values
    int k = 5; // k-mer length
    int b = 4; // number of bins/buckets to assign canonical k-mers to

    // parse arguments
    if (argc < 2) {
        fprintf(stderr, "arguments: fasta file, k (default 5, <=32), bin number (default 4)");
        exit(1);
    }

    if (argc > 2) {
        k = atoi(argv[2]);
        if (0 == k || k > 32) {
            fprintf(stderr, "ERROR: k must be in [1,32].\n");
            exit(1);
        }
    }

    if (argc > 3) {
        b = atoi(argv[3]);
    }

    // initialize array
    int bins[b];
    for (int i = 0; i < b; i++) {
        bins[i] = 0;
    }

    /* process FASTA file */
    FASTAFILE* ffp;
    char* seq;
    char* name;
    int L;
    ffp = OpenFASTA(argv[1]);
    while (ReadFASTA(ffp, &seq, &name, &L)) {

        if (
            //***
            //*** distribute canonical k-mers to bins
            //***

            process_string(seq, k, bins, b)
            // process_string_std(seq,k,bins,b)
        ) {
            exit(1);
        }

        free(seq);
        free(name);
    }
    CloseFASTA(ffp);

    // output distribution to bins
    int sum = 0;
    for (int i = 0; i < b; i++) {
        printf("%d\n", bins[i]);
        sum += bins[i];
    }
    // printf("\nSUM: %d\n",sum);

    return 0;
}
