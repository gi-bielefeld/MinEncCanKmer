#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fasta.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdbool.h>


// print a k-mer in binary format
void printkmer(u_int64_t u){
		// print kmer
 		u_int64_t t = pow(2, 63);   // t is the max number that can be represented
	
		for(t; t>0; t = t/2){       // t iterates through powers of 2
			if(u >= t){             // check if u can be represented by current value of t
				u -= t;
				printf("1");        // if so, add a 1
			}
			else {
				printf("0");        // if not, add a 0
			}
		}
		
		printf("\n");	
}

// extract all k-mers, encode by standard 2-bit encoding, assign to b bins
int process_string_std(char* s, int k, int* bins, int b){

	u_int64_t kmer=0;
	u_int64_t rckmer=0;
	u_int64_t rem = 3;
	rem<<= (2*k);
	u_int64_t pos1 = 1;
	pos1<<= (2*k-2);
	u_int64_t pos2 = 1;
	pos2<<= (2*k-1);
	u_int64_t maxrank = powl(2,2*k);
	char warned=0; //warn first time, an unsupported character is skipped

	for (int i=0;i<strlen(s);i++){

		//update kmer by current character
		kmer=kmer<<2;
		switch(s[i]) {
			case 'A': break;
			case 'C': kmer+=1; break;
			case 'G': kmer+=2; break;
			case 'T': kmer+=3; break;
			default: if(!warned) {fprintf(stderr, "Warning: unsupported character replaced by A: %c\n", s[i]); warned=1; } break;
		}

		// clear unused bits
 		kmer &= ~rem;
		
		// update reverse complement kmer
		rckmer=rckmer>>2;
		switch(s[i]) {
			case 'T': break;
			case 'G': rckmer|=pos1; break;
			case 'C': rckmer|=pos2; break;
			case 'A': rckmer|=pos1|pos2; break;
			default: rckmer|=pos1|pos2; break;
		}
		
	   if (i>=k-1){
		   
			
			u_int64_t can = kmer < rckmer ? kmer : rckmer;
			
			// assign to bin
			u_int64_t p = b * can / maxrank ;
// 			int p = can % t ;
 			bins[p]++;
	   }

	}
	
	return 0;
}


// extract all k-mers, encode by enc^r_c, assign to b bins
int process_string(char* s, int k, int* bins, int b){


	
	char warned=0; //warn first time, an unsupported character is skipped
	u_int64_t kmer=0;
	u_int64_t rckmer=0;

	int max=8*sizeof(kmer);
	int offset = max-2*k;

	// precompute masks
	u_int64_t allones = powl(2,max)-1;

    u_int64_t posmasks[max+1];
	posmasks[max]=1;
	for(int i=max-1; i>=0; i-- ){
		posmasks[i]=posmasks[i+1]<<1;
	}

	// 1 .. 1 0 .. 0
	u_int64_t onemasks[max+1];
	onemasks[max]=allones;
	for(int i=max-1; i>=0; i-- ){
		onemasks[i]=(onemasks[i+1]<<1);
	}

	// 0 .. 0 1 .. 1
	u_int64_t zeromasks[max+1];
	zeromasks[0]=allones;
	for(int i=1; i<=max; i++ ){
		zeromasks[i]=zeromasks[i-1]>>1;
	}

	// 0 .. 0 1 .. 1 0 .. 0
	u_int64_t remaindermasks[k+1];
	remaindermasks[0]=allones;
	for(int i=1; i<=k; i++ ){
		remaindermasks[i]=zeromasks[i+offset]&onemasks[max-i];
	}
	
	u_int64_t maxrank = powl(2,2*k-1);
	
	if(k%2==0){
		//even k -> need to consider palindromes
		maxrank=powl(4,k)/2;
		// this is not the actual max rank. It is chosen such that the bucket sizes are a bit smaller -- such that palindromic k-mers (which do not appear in pairs) contribute as much as the other canonical k-mers
	}
	

	// K:
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
	char replace1[16] = {0,0,0,0,1,0,0,0,1,0,0,0,0,1,1,0};
	char replace2[16] = {1,1,1,0,0,1,0,1,0,0,1,1,0,0,0,1};
	char replace3[16] = {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
	char replace4[16] = {0,1,0,0,0,1,0,0,1,0,1,1,0,1,0,0};
	char reverse[16] =  {0,0,0,0,0,0,0,1,0,0,1,1,0,1,1,1};

	// compute encoding where only setting the bits accorodung to specifying case and subtracting gaps is missing
	u_int64_t initialhashed(u_int64_t *kmer, int offset, int l){

		// pick a precomputed mask consisting of l trailing 1s and 0s else
		// do and AND of that mask and the original k-mer to get the new right part
		u_int64_t right = (*kmer) & zeromasks[offset+2*k-l]; //l trailing ones
		
		// complement/invert the right part by xor with zeromask
		// /!\ TO SAVE COMPUTATION TIME, THIS STEP CAN BE SKIPPED. The resulting encodung will be different bbut still a coorect minimal encoding.
		right = right ^ zeromasks[offset+2*k-l];

		// no remainder left?
		if((l+2)>=k){ return(right); }
		
		// pick a precomputed mask consisting of ones in the middle to get the remainder
		u_int64_t remainder = (*kmer) & remaindermasks[l+2];
		
		// shift remainder two bits to the right
		remainder=remainder>>2;
			
		// do OR of left and middle part
		return (remainder | right);
	}

	// compute enc_r_c
	u_int64_t hash(u_int64_t *kmer, u_int64_t *rckmer){
		
		// hash
		u_int64_t kmerhash;
		
		// get length of symmetric pre/suffix
		u_int64_t sym = (*kmer) ^ (*rckmer);
		int l = __builtin_ctzll(sym)/2*2;


		if (l<k-1){ // not just single character in the middle
			
			// get the first two asymmetric characters, i.e. 2x2 bits
			char pattern=0;
			if((*kmer)&posmasks[offset+l+1]){ pattern+=8;}
			if((*kmer)&posmasks[offset+l+2]){ pattern+=4;}
			if((*kmer)&posmasks[offset+2*k-l-1]){ pattern+=2;}
			if((*kmer)&posmasks[offset+2*k-l]){ pattern+=1;}
			
			if(reverse[pattern]){
				kmerhash = initialhashed(rckmer,offset,l);
			} else {
				kmerhash = initialhashed(kmer,offset,l);
			}

			// set positions l+1, l+2, l+3 and l+4 according to *-pair-encoding
			if(replace1[pattern]){ kmerhash|=posmasks[offset+l+1];}
			if(replace2[pattern]){ kmerhash|=posmasks[offset+l+2];}
			if(replace3[pattern]){ kmerhash|=posmasks[offset+l+3];}
			if(replace4[pattern]){ kmerhash|=posmasks[offset+l+4];}

			
		} else if (l>=k) { // palindrome -> nothing to do
			l=k;
			kmerhash = initialhashed(kmer,offset,l);
		} else { // single character in the middle
			kmerhash = initialhashed(kmer,offset,l);
			// set the bits accordingly
			// A=00 -> 0
			// C=01 -> 1
			// G=10 -> 1
			// T=11 -> 0
			if((*kmer)&posmasks[offset+l+1] && !(*kmer)&posmasks[offset+l+2]){ kmerhash|=posmasks[offset+l+2];} //rc
			if((*kmer)&posmasks[offset+l+2]){ kmerhash|=posmasks[offset+l+2];}
		} 
		
		
		// subtract gaps
		// 2*(k//2-l-1) ones followed by k-2 zeros
		if(l<=k-4){
			u_int64_t gaps=zeromasks[max-(2*(k/2-l/2-1))];
			gaps=gaps<<(2*((k+1)/2)-1);
			kmerhash-=gaps;
		}

		// subtract gap in code due to specifying middle position
		if (k%2==1 && kmerhash>=pow(4,(k/2+1))){
			kmerhash -= 2*pow(4,k/2);
		}
		
		return kmerhash;
	}

	// go through the string, process each k-mer by shift/update, enc_r_c, and assign to bins
	for (int i=0;i<strlen(s);i++){
		

		//update kmer by current character
		kmer=kmer<<2;
		switch(s[i]) {
			case 'A': break;
			case 'C': kmer+=1; break;
			case 'G': kmer+=2; break;
			case 'T': kmer+=3; break;
			default: if(!warned) {fprintf(stderr, "Warning: unsupported character replaced by A: %c\n", s[i]); warned=1; } break;
		}
		
		// clear unused bits
 		kmer &= ~posmasks[offset];
 		kmer &= ~posmasks[offset-1];

		
		// update reverse complement kmer
		rckmer=rckmer>>2;
		switch(s[i]) {
			case 'T': break;
			case 'G': rckmer|=posmasks[offset+2]; break;
			case 'C': rckmer|=posmasks[offset+1]; break;
			case 'A': rckmer|=posmasks[offset+2]|posmasks[offset+1]; break;
			default:  rckmer|=posmasks[offset+2]|posmasks[offset+1]; break;
		}
		

	   if (i>=k-1){
		   
	   
			// hash
		    u_int64_t kmerhash = hash(&kmer, &rckmer);
			
			// assign to bin
			u_int64_t p = b*kmerhash/maxrank;
			
			if (p==b) {
				// palindromes comes not in pairs like the other canonical k-mers. In order to ensure a equal distribution, bucket sizes are chosen as if there were only half as many palindromes possible. This results in ranks that would correspond to a max+1st bin. These k-mers are assigned to bin 0 where the palindromes are "missing".
				bins[0]++;
			} else {
				bins[p]++;
			}
				
			
		}

	}
	
	return 0;
}


int main(int argc, char* argv[]) {
   
	// default values
	int k=5; // k-mer length
	int b=4; // number of bins/buckets to assign canonical k-mers to

	// parse arguments
	if (argc<2){
		fprintf(stderr, "arguments: fasta file, k (default 5, smaller 32), bin number (default 4)");
		exit(1);
	}
	
	if (argc>2){
		k=atoi(argv[2]);
		if(k>31){
			fprintf(stderr, "ERROR: k must be smaller than 32.\n");
			exit(1);
		}
	}
	
	if (argc>3){
		b=atoi(argv[3]);
	}


	
	// initialize array
	int bins[b];
	for (int i=0;i<b;i++){
		bins[i]=0;
	}
	
	/* process FASTA file */
	FASTAFILE *ffp;
	char *seq;
	char *name;
	int   L;
	ffp = OpenFASTA(argv[1]);
	while (ReadFASTA(ffp, &seq, &name, &L)) {

		
		if (
		
			//***
			//*** distribute canonical k-mers to bins
			//***
			
			process_string(seq,k,bins,b)
// 			process_string_std(seq,k,bins,b)
			
		){exit(1);}

		free(seq);
		free(name);
	}
	CloseFASTA(ffp);
  
	
	// output distribution to bins
	int sum=0;
	for (int i=0;i<b;i++){
 		printf("%d\n",bins[i]);
		sum+=bins[i];
	}
//  	printf("\nSUM: %d\n",sum);
   
	return 0;
   
}
