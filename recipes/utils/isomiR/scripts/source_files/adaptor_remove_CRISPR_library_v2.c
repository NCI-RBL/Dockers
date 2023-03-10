
// 2015 Jan 6th  Rewrite the old version to analyze data from illumina small RNA kit.

// search 3' adaptor sequnce and remove it.

// based on the percentage of match- for example 0.8 means 1 mismatch is allowed out of 5nt

// find the first one then remove the rest.

// maybe I need to add a feature to remove all contaminat

// maybe I need to give stat about length distribution

// 2016 Jan Find out that the first 5nt of hsa-miR-1 is excact the same as adaptor. Changed the pass-rate to 0.9 and minimal length to 10. 

// 2016 Jan 20th, change the minimal length of reads to 10nt -> to analyze the trimmed reads of miR-27a/b/s

// 2016 May 10th, change the adaptor sequence to NEB, keep the min length as 10nt

// 2016 Nov 14th. Completely rewrite for Lisheng to analyze the deep seq result of CRISPR screen. 1. find 5' adaptor. 2. take 20nt sgRNA seq. 3. check if the following sequence are tracrRNA.

// 2016 Nov 18th. output to file. also output reads w/o adaptors.



// v2 2016 Nov 19th. Change the mechanism to search for adaptors. 1. search for 5' adaptor.  2. then search for 3' adaptor. 3. determine the length of inserts. 4. trim the reads to 20nt long.

// filter away reads shorter than 17nt!  reads with a size between 17nt and 23nt will be analyzed! (17 <= n <= 23)



#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define Bool	int
#define	TRUE	1
#define	FALSE	0

#define	LEN_NAME	100
#define	LEN_SEQ		100
#define	LEN_TAG		15

#define ADAPTOR_U6              "GTGGAAAGGACGAAACACCG"          //  Of note: the last nt is G for U6 transcription start!
#define ADAPTOR_tracrRNA		"GTTTTAGAGCTAGAAATAGC"
#define	MISMATCH_5          1
#define MISMATCH_3          1

	
struct fastq 
{
	char 	name[LEN_NAME+1];
	char 	sequence[LEN_SEQ+1];
	char	tag[LEN_TAG+1];
	char	score[LEN_SEQ+1];
};

char*	read_fastq(struct fastq*);

void	print_fastq(const struct fastq*, FILE*);

char*	found_adaptors(char*, char*, int);	/* return TURE if adaptor is found, otherwise return FALSE */

Bool	adaptor_remove( struct fastq* );  // removed return TRUE, otherwise return FALSE

int strcmp_mis ( char* s1, char* s2, int n);

long int 	num_read = 0,  contain_5adaptor = 0, contain_both_adaptors=0, right_length=0, shorter=0, longer=0, too_long=0, too_short=0, num_output=0;

int main(int argc, char* argv[])
{
	
    FILE                *fp_output, *fp_no_adaptor;
    
    if (argc != 3)
	{
		printf("\nadaptor_remove.exe output.fastq no_adaptor.fastq < input.fastq\n");
		return 0;
	}
    
    fp_output = fopen(argv[1], "w");
    fp_no_adaptor = fopen(argv[2], "w");

	struct fastq	read_temp;
	struct fastq 	*p_read = &read_temp;	/* p_read point to the current read */

	fprintf(stderr, "\n U6 Adaptor sequence is %s, allowing %d mismatch(es). TracrRNA Adaptor sequence is %s, allowing %d mismatch(es)\n", ADAPTOR_U6, MISMATCH_5, ADAPTOR_tracrRNA, MISMATCH_3) ;
	
	while (read_fastq(p_read) != NULL )
	{	
		num_read++;
		
		if ( adaptor_remove(p_read) )   {print_fastq(p_read, fp_output); num_output++;}
		
        else                            print_fastq(p_read, fp_no_adaptor);
    }
	
	
	fprintf(stderr, "Total number of reads:			 %ld\n", num_read);
	fprintf(stderr, "number of reads have U6_adaptor :       %ld, that is %ld%% of total reads\n", contain_5adaptor , contain_5adaptor*100/num_read);
    fprintf(stderr, "Among these, %ld reads have Tracer_adaptor, that is %ld%% of last step and %ld%% of total reads\n", contain_both_adaptors, contain_both_adaptors*100/contain_5adaptor , contain_both_adaptors*100/num_read);
    fprintf(stderr, "%ld reads are too short (<17nt), %ld reads are too long (>23nt). These reads are discarded.\n", too_short, too_long);
    fprintf(stderr, "In the rest of reads, %ld reads are 20nt, %ld reads are shorter than 20nt, %ld reads are longer than 20nt\n", right_length, shorter, longer);
    fprintf(stderr, "Final output are %15ld reads. That is %ld%% of reads with both adaptors and %ld%% of total reads\n", num_output, num_output*100/contain_both_adaptors, num_output*100/num_read);
    
    
	fclose(fp_output);
    fclose(fp_no_adaptor);
    
    return 0;
}

char* read_fastq(struct fastq* read)
{
	if (gets(read->name) == NULL) return NULL;
	gets(read->sequence);
	gets(read->tag);
	gets(read->score);
	
	return read->name;
}

void print_fastq(const struct fastq* read, FILE* fp)
{
	fprintf(fp, "%s\n", read->name);
	fprintf(fp,"%s\n", read->sequence);
	fprintf(fp,"%s\n", read->tag);
	fprintf(fp,"%s\n", read->score);
}


		
Bool	adaptor_remove( struct fastq* read)
{
	
	char*   pos1 = NULL;
    char*   pos2 = NULL;
    int     length_5adaptor;
    int     length_sgRNA;
    char    sgRNA_temp[21];           // To hold the sgRNA sequence.
    
    
    
    if ((pos1 = found_adaptors(read->sequence, ADAPTOR_U6, MISMATCH_5)) == NULL) return FALSE;
    else
    {
        contain_5adaptor++;
        length_5adaptor = strlen(ADAPTOR_U6);
        pos1 += length_5adaptor;                    //pos1 point to the nt after 5_adaptor
        
        if ((pos2 = found_adaptors(pos1, ADAPTOR_tracrRNA, MISMATCH_3)) == NULL) return FALSE;
        contain_both_adaptors++;                    // pos2 point to the first nt of 3_adaptor
        
        length_sgRNA = pos2 - pos1;
        if ( length_sgRNA > 23) { too_long++; return FALSE;}
        if ( length_sgRNA < 17) { too_short++; return FALSE;}
        if ( length_sgRNA == 20) right_length++;
        else if (length_sgRNA < 20) shorter++;
        else { longer++; pos2 = pos1 + 20; length_sgRNA = 20;}             // trim long reads to 20nt long
            
        *pos2 = '\0';
        strcpy (sgRNA_temp, pos1);
        strcpy (read->sequence, sgRNA_temp);
        
         *(read->score + length_sgRNA) = '\0';  // trim the score to the right length
    
        return TRUE;
    }
    
}

    
int strcmp_mis ( char* s1, char* s2, int n)     /* s1 and s2 are the pointers to the strings to compare, n is the numbr of mismatches allowed */
                                                /* for return, -1 - no match, mismatches are more than n,  m>=0 - remaining mismatches allowrance */
                                                /* s1 and s2 both need a \0 at the the ends */
{
    while ( n >= 0 )
    {
        if ( *s1 - *s2 )    n--;
        if ( *s1 )          s1++;
        if ( *s2 )          s2++;
        if (*s2 == '\0') return n;             // search finish when the end of key is researched
    }
    return n;
}

char*   found_adaptors (char* s1, char* s2, int n )     /* Returns a pointer to the first occurrence of str2 in str1, or a null pointer if str2 is not part of str1. */
                                                    /* s1 is the string to be searched and s2 is the key. */
                                                    /* n is the mismatch allowrance */
                                                    /* strstr_mis search the whole s1, and report the FIRST matching of s2 within the lowest mismatch number */

{
    char * pos = NULL;
    int min_mis = -1, mismatch;
    
    while (*s1)
    {
        mismatch = strcmp_mis(s1, s2, n);
        if (mismatch > min_mis) {min_mis = mismatch; pos = s1;}
        s1++;
    }
    
    return pos;
}






