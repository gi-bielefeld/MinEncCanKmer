#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fasta.h>
#include <stdlib.h>
#include <ctype.h>



int process_string_std(char* s, int k, int* threads, int t){
  


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
			
			// assign to thread
			u_int64_t p = t * can / maxrank ;
// 			int p = can % t ;
 			threads[p]++;
	   }

	}
	
	return 0;
}



int process_string(char* s, int k, int* threads, int t){


	
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

	u_int64_t onemasks[max+1];
	onemasks[max]=allones;
	for(int i=max-1; i>=0; i-- ){
		onemasks[i]=(onemasks[i+1]<<1);
		onemasks[i]=onemasks[i]&~(posmasks[offset+1]);
		onemasks[i]=onemasks[i]|(posmasks[offset]);
	}


	u_int64_t zeromasks[max+1];
	zeromasks[0]=allones;
	for(int i=1; i<=max; i++ ){
		zeromasks[i]=zeromasks[i-1]>>1;
	}

	
	u_int64_t maxrank = powl(2,2*k-1);
	
	if(k%2==0){
		//even k -> need to consider palindromes
		maxrank=powl(4,k)/2;
		// this is not the actual max rank. It is chosen such that the bucket sizes are a bit smaller -- such that palindromic k-mers (which do not appear in pairs) contribute as much as the other canonical k-mers
	}
	
	// * 0 A..A -> 00..0
	// * 1 A..C -> 00..1
	// * 2 A..G -> 01..0
	// # 3 palindrome A..T -> 11..1
	// * 4 C..A -> 01..1
	// * 5 C..C -> 10..0
	// # 6 palindrome C..G -> 11..0
	// # 7 C..T -> A..G -> 01..0
	// * 8 G..A -> 10..1
	// # 9 palindrome G..C -> 00..1
	// # 10 G..G -> C..C -> 10..0
	// # 11 G..T -> A..C -> 00..1
	// # 12 palindrome T..A -> 00..0
	// # 13 T..C -> G..A -> 10..1
	// # 14 T..G -> C..A -> 01..1
	// # 15 T..T -> A..A -> 00..0
	char replace1[16] = {0,1,0,1,1,0,0,0,1,1,0,1,0,1,1,0};
	char replace2[16] = {0,0,1,1,1,0,1,1,0,0,0,0,0,0,1,0};
	char replace3[16] = {0,0,0,1,0,1,1,0,1,0,1,0,0,1,0,0};
	char reverse[16] =  {0,0,0,0,0,0,0,1,0,0,1,1,0,1,1,1};


	
u_int64_t initialshortened(u_int64_t *kmer, int offset, int l){

		const u_int64_t shiftedkmer = (*kmer)>>1;
	
		// pick a precomputed mask consisting of l trailing 1s and 0s else
		// do and AND of that mask and the original k-mer to get the new right part
		u_int64_t right = (*kmer) & zeromasks[offset+2*k-l]; //l trailing ones

		// pick a precomputed mask consisting of a 0 , 1s, and l trailing 0s
		// do and AND of that mask and the right-shifted k-mer to get the left part
		u_int64_t left = shiftedkmer & onemasks[offset+2*k-l]; //l trailing zeros

		// do and OR of  left (unshifted) and right (shifted) part
		return (left | right);
}


u_int64_t shorten(u_int64_t *kmer, u_int64_t *rckmer){
		   
			
			// shorten
		   u_int64_t shortenedkmer;
		   
			// get length of symmetric pre/suffix
			u_int64_t sym = (*kmer) ^ (*rckmer);
			int l = __builtin_ctzll(sym)/2*2;
			
			// l<k-1 -> switch to inner encoding
			// l==k-1 -> single charachter in the middle
			// l==64 -> palindrome (set l to k-2!?)
			
			char pattern=0;
			
			//palindrome. Set l such that coordinates match with usual "switch to inner" mode
			if (l>k-1){l=k-2;} 

			if (l<k-1){ // not just single character in the middle
				
				// get the first two asymmetric characters, i.e. 2x2 bits
				if((*kmer)&posmasks[offset+l+1]){ pattern+=8;}
				if((*kmer)&posmasks[offset+l+2]){ pattern+=4;}
				if((*kmer)&posmasks[offset+2*k-l-1]){ pattern+=2;}
				if((*kmer)&posmasks[offset+2*k-l]){ pattern+=1;}
				
				if(reverse[pattern]){
					shortenedkmer = initialshortened(rckmer,offset,l);
				} else {
					shortenedkmer = initialshortened(kmer,offset,l);
				}

				// set positions l+2, l+3 and 2k-l according to *-pair-encoding
				//setpos(&shortenedkmer,offsetl,offsetml,pattern);
				if(replace1[pattern]){
					shortenedkmer|=posmasks[offset+2*k-l];
				} else {
					shortenedkmer&=~posmasks[offset+2*k-l];
				}
				if(replace2[pattern]){
					shortenedkmer|=posmasks[offset+l+3];
				} else {
					shortenedkmer&=~posmasks[offset+l+3];
				}
				if(replace3[pattern]){
					shortenedkmer|=posmasks[offset+l+2];
				} else {
					shortenedkmer&=~posmasks[offset+l+2];
				}

			} else { // single character in the middle
				
				shortenedkmer = initialshortened(kmer,offset,l);
			
				// set the single bits according to encoding of *-case
				//A,G -> do nothing
				//C->1
				if(!((*kmer)&posmasks[offset+l+1]) && (*kmer)&posmasks[offset+l+2]){
					shortenedkmer|=posmasks[offset+2*k-l];
				} else{
					//T->0
					if((*kmer)&posmasks[offset+l+1] && (*kmer)&posmasks[offset+l+2]){
						shortenedkmer&=~posmasks[offset+2*k-l];
					}
				}
				
			} 

			// pick a precomputed mask consisting of l+1 1s and then 0s
			// OR with first mask to set 1s for left outer characters
			shortenedkmer|=onemasks[offset+l+1];
			shortenedkmer&=~onemasks[offset];

			
			// palindrome
			if(pattern==9 || pattern==12){
				// set first half to 0
				shortenedkmer&=zeromasks[offset+k];
				//  set very first bit to 1 to indicate palindrome
				shortenedkmer|=posmasks[offset+1];
				// for other palindromes (patterns 3 and 6), nothing has to be done
			}

			return shortenedkmer;
}


	
	
	
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
		   
			// shorten
		    u_int64_t shortenedkmer = shorten(&kmer, &rckmer);
			
			// assign to thread
			u_int64_t p = t*shortenedkmer/maxrank;

			if (p==t) {
				// palindromes comes not in pairs like the other canonical k-mers. In order to ensure a equal distribution, bucket sizes are chosen as if there were only half as many palindromes possible and the other palindromes are assigned to the last bucket
				threads[p-1]++;
			} else {
				threads[p]++;
			}
				
//    		printf("%llu\n",p);
//  		printf("%llu\n",shortenedkmer);
			
		}

	}
	
	return 0;
}


int main(int argc, char* argv[]) {
   
	// default values
	int k=5; // k-mer length
	int t=4; // number of threads/buckets to assign canonical k-mers to

	// parse arguments
	if (argc<2){
		fprintf(stderr, "arguments: fasta file, k (default 5, smaller 32), t (default 4)");
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
		t=atoi(argv[3]);
	}


	
	// initialize array
	int threads[t];
	for (int i=0;i<t;i++){
		threads[i]=0;
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
			//*** distribute canonical k-mers to threads
			//***
			
			process_string(seq,k,threads,t)
// 			process_string_std(seq,k,threads,t)
			
		){exit(1);}

		free(seq);
		free(name);
	}
	CloseFASTA(ffp);
  

//  	exit(0);
	
	// output distribution to threads
	int sum=0;
	for (int i=0;i<t;i++){
 		printf("%d ",threads[i]);
		sum+=threads[i];
	}
 	printf("\nSUM: %d\n",sum);
   
	return 0;
   
}
