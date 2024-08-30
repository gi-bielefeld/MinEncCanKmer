#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fasta.h"

// =================================================================================================
//     Helper Functions
// =================================================================================================

// print a k-mer in binary format
void printkmer(u_int64_t u)
{
    // print kmer
    u_int64_t t = pow(2, 63); // t is the max number that can be represented

    for (t; t > 0; t = t / 2) { // t iterates through powers of 2
        if (u >= t) { // check if u can be represented by current value of t
            u -= t;
            printf("1"); // if so, add a 1
        } else {
            printf("0"); // if not, add a 0
        }
    }

    printf("\n");
}

// =================================================================================================
//     Minimal Canonical Encoding Functions
// =================================================================================================

typedef struct {
    int k;

    u_int64_t allones;
    int max;
    int offset;

    // max + 1, since max = 8 * sizeof(u_int64_t)
    u_int64_t posmasks[65];
    u_int64_t onemasks[65];
    u_int64_t zeromasks[65];
    u_int64_t remaindermasks[34]; // max k + 2
} Bitmasks;

// Fixed bitmasks to detect the pattern of a specifying pair.
// Because this is C, it is a bit tricky to make them constant and in the right scope,
// so for simplicity, we just make them global.
// R:
// * 0 A..A -> 0110
// * 1 A..C -> 0101
// * 2 A..G -> 0100
// # 3 palindrome A..T
// * 4 C..A -> 1000
// * 5 C..C -> 0111
// # 6 palindrome C..G
// # 7 C..T -> A..G -> 0100
// * 8 G..A -> 1001
// # 9 palindrome G..C
// # 10 G..G -> C..C -> 0111
// # 11 G..T -> A..C -> 0101
// # 12 palindrome T..A
// # 13 T..C -> G..A -> 1001
// # 14 T..G -> C..A -> 1000
// # 15 T..T -> A..A -> 0110
static const char replace1[16] = { 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0 };
static const char replace2[16] = { 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1 };
static const char replace3[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };
static const char replace4[16] = { 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0 };
static const char reverse[16] = { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1 };

void initialize_bitmasks(Bitmasks* bm, int k)
{
    bm->k = k;
    bm->max = 8 * sizeof(u_int64_t);
    bm->offset = bm->max - 2 * k;

    // Precompute masks
    bm->allones = powl(2, bm->max) - 1;

    // Initialize posmasks
    bm->posmasks[bm->max] = 1;
    for (int i = bm->max - 1; i >= 0; i--) {
        bm->posmasks[i] = bm->posmasks[i + 1] << 1;
    }

    // Initialize onemasks
    // 1 .. 1 0 .. 0
    bm->onemasks[bm->max] = bm->allones;
    for (int i = bm->max - 1; i >= 0; i--) {
        bm->onemasks[i] = (bm->onemasks[i + 1] << 1);
    }

    // Initialize zeromasks
    // 0 .. 0 1 .. 1
    bm->zeromasks[0] = bm->allones;
    for (int i = 1; i <= bm->max; i++) {
        bm->zeromasks[i] = bm->zeromasks[i - 1] >> 1;
    }

    // Initialize remaindermasks
    bm->remaindermasks[0] = bm->allones;
    for (int i = 1; i <= k; i++) {
        bm->remaindermasks[i] = bm->zeromasks[i + bm->offset] & bm->onemasks[bm->max - i];
    }
}

// compute encoding where only setting the bits accorodung to specifying case and subtracting gaps is missing
u_int64_t encode_prime(Bitmasks const* bm, u_int64_t kmer, int offset, int l)
{
    int k = bm->k;

    // pick a precomputed mask consisting of l trailing 1s and 0s else
    // do and AND of that mask and the original k-mer to get the new right part
    u_int64_t right = kmer & bm->zeromasks[offset + 2 * k - l]; //l trailing ones

    // complement/invert the right part by xor with zeromask
    // /!\ TO SAVE COMPUTATION TIME, THIS STEP CAN BE SKIPPED. The resulting encodung will be different bbut still a coorect minimal encoding.
    right = right ^ bm->zeromasks[offset + 2 * k - l];

    // no remainder left?
    if ((l + 2) >= k) {
        return (right);
    }

    // pick a precomputed mask consisting of ones in the middle to get the remainder
    u_int64_t remainder = kmer & bm->remaindermasks[l + 2];

    // shift remainder two bits to the right
    remainder = remainder >> 2;

    // do OR of left and middle part
    return (remainder | right);
}

// compute enc_r_c
u_int64_t encode(Bitmasks const* bm, u_int64_t kmer, u_int64_t rckmer)
{
    // hash
    u_int64_t kmerhash;
    int k = bm->k;

    // get length of symmetric pre/suffix
    u_int64_t sym = kmer ^ rckmer;
    int l = __builtin_ctzll(sym) / 2 * 2;

    if (l < k - 1) { // not just single character in the middle

        // get the first two asymmetric characters, i.e. 2x2 bits
        char pattern = 0;
        if (kmer & bm->posmasks[bm->offset + l + 1]) {
            pattern += 8;
        }
        if (kmer & bm->posmasks[bm->offset + l + 2]) {
            pattern += 4;
        }
        if (kmer & bm->posmasks[bm->offset + 2 * k - l - 1]) {
            pattern += 2;
        }
        if (kmer & bm->posmasks[bm->offset + 2 * k - l]) {
            pattern += 1;
        }

        if (reverse[pattern]) {
            kmerhash = encode_prime(bm, rckmer, bm->offset, l);
        } else {
            kmerhash = encode_prime(bm, kmer, bm->offset, l);
        }

        // set positions l+1, l+2, l+3 and l+4 according to *-pair-encoding
        if (replace1[pattern]) {
            kmerhash |= bm->posmasks[bm->offset + l + 1];
        }
        if (replace2[pattern]) {
            kmerhash |= bm->posmasks[bm->offset + l + 2];
        }
        if (replace3[pattern]) {
            kmerhash |= bm->posmasks[bm->offset + l + 3];
        }
        if (replace4[pattern]) {
            kmerhash |= bm->posmasks[bm->offset + l + 4];
        }

    } else if (l >= k) { // palindrome -> nothing to do
        l = k;
        kmerhash = encode_prime(bm, kmer, bm->offset, l);
    } else { // single character in the middle
        kmerhash = encode_prime(bm, kmer, bm->offset, l);
        // set the bits accordingly
        // A=00 -> 0
        // C=01 -> 1
        // G=10 -> 1
        // T=11 -> 0
        if ((kmer & bm->posmasks[bm->offset + l + 1]) && !(kmer & bm->posmasks[bm->offset + l + 2])) {
            kmerhash |= bm->posmasks[bm->offset + l + 2];
        } //rc
        if (!(kmer & bm->posmasks[bm->offset + l + 1]) && (kmer & bm->posmasks[bm->offset + l + 2])) {
            kmerhash |= bm->posmasks[bm->offset + l + 2];
        }
    }

    // subtract gaps
    // 2*(k//2-l-1) ones followed by k-2 zeros
    if (l <= k - 4) {
        u_int64_t gaps = bm->zeromasks[bm->max - (2 * (k / 2 - l / 2 - 1))];
        gaps = gaps << (2 * ((k + 1) / 2) - 1);
        kmerhash -= gaps;
    }

    // subtract gap in code due to specifying middle position
    if (k % 2 == 1 && kmerhash >= pow(4, (k / 2 + 1))) {
        kmerhash -= 2 * pow(4, k / 2);
    }

    return kmerhash;
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

    u_int64_t maxrank = powl(2, 2 * k - 1);
    if (k % 2 == 0) {
        //even k -> need to consider palindromes
        maxrank = powl(4, k) / 2;
        // this is not the actual max rank. It is chosen such that the bucket sizes are a bit smaller -- such that palindromic k-mers (which do not appear in pairs) contribute as much as the other canonical k-mers
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
        kmer &= ~bm.posmasks[bm.offset];
        kmer &= ~bm.posmasks[bm.offset - 1];

        // update reverse complement kmer
        rckmer = rckmer >> 2;
        switch (s[i]) {
        case 'T':
            break;
        case 'G':
            rckmer |= bm.posmasks[bm.offset + 2];
            break;
        case 'C':
            rckmer |= bm.posmasks[bm.offset + 1];
            break;
        case 'A':
            rckmer |= bm.posmasks[bm.offset + 2] | bm.posmasks[bm.offset + 1];
            break;
        default:
            rckmer |= bm.posmasks[bm.offset + 2] | bm.posmasks[bm.offset + 1];
            break;
        }

        // Not yet done filling the first kmer
        if (i < k - 1) {
            continue;
        }

        // hash
        u_int64_t kmerhash = encode(&bm, kmer, rckmer);

        // assign to bin
        u_int64_t p = b * kmerhash / maxrank;

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
    u_int64_t rem = 3;
    rem <<= (2 * k);
    u_int64_t pos1 = 1;
    pos1 <<= (2 * k - 2);
    u_int64_t pos2 = 1;
    pos2 <<= (2 * k - 1);
    u_int64_t maxrank = powl(2, 2 * k);
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
            // 			int p = can % t ;
            bins[p]++;
        }
    }

    return 0;
}

// =================================================================================================
//     Main
// =================================================================================================

int main(int argc, char* argv[])
{

    // default values
    int k = 5; // k-mer length
    int b = 4; // number of bins/buckets to assign canonical k-mers to

    // parse arguments
    if (argc < 2) {
        fprintf(stderr, "arguments: fasta file, k (default 5, smaller 32), bin number (default 4)");
        exit(1);
    }

    if (argc > 2) {
        k = atoi(argv[2]);
        if (k > 31) {
            fprintf(stderr, "ERROR: k must be smaller than 32.\n");
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
            // 			process_string_std(seq,k,bins,b)

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
    //  	printf("\nSUM: %d\n",sum);

    return 0;
}
