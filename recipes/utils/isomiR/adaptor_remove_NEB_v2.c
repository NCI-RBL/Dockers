
// 2015 Jan 6th  Rewrite the old version to analyze data from illumina small RNA kit.

// search 3' adaptor sequnce and remove it.

// based on the percentage of match- for example 0.8 means 1 mismatch is allowed out of 5nt

// find the first one then remove the rest.

// maybe I need to add a feature to remove all contaminat

// maybe I need to give stat about length distribution

// 2016 Jan Find out that the first 5nt of hsa-miR-1 is excact the same as adaptor. Changed the pass-rate to 0.9 and minimal length to 10. 

// 2016 Jan 20th, change the minimal length of reads to 10nt -> to analyze the trimmed reads of miR-27a/b/s

// 2016 May 10th, change the adaptor sequence to NEB, keep the min length as 10nt

// 2017 Dec 4th, re-compile the current version to v1.exe

// 2023 Mar 22nd, allowing users to input min length; change the name to v2.

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define Bool	int
#define	TRUE	1
#define	FALSE	0

#define	LEN_NAME	150
#define	LEN_SEQ		150
#define	LEN_TAG		150

#define ADAPTOR		"AGATCGGAAGAGCACACGTCT"
#define	PASS_RATE           0.9
#define MIN_MATCH           10

	
struct fastq 
{
	char 	name[LEN_NAME+1];
	char 	sequence[LEN_SEQ+1];
	char	tag[LEN_TAG+1];
	char	score[LEN_SEQ+1];
};

char*	read_fastq(struct fastq*);

void	print_fastq(const struct fastq*);

Bool	found_adaptors(char*);	/* return TURE if adaptor is found, otherwise return FALSE */

Bool	adaptor_remove( struct fastq* );  // removed return TRUE, otherwise return FALSE

int main(int argc, char* argv[])		
{
	if (argc != 2) 
	{
		printf("\n adaptor_remove_NEB_v2.exe min_length_pass_filter  < input.fastq > output.fastq\n"); 
		return 0;
	}
	
	int         min_length_pass_filter = atoi (argv[1]); // new code added - march 2023

    long int      num_read = 0,  contain_adaptors = 0, too_short = 0;

	struct fastq	read_temp;
	struct fastq 	*p_read = &read_temp;	/* p_read point to the current read */

	fprintf(stderr, "\n3' Adaptor sequence is %s and passing_rate is %f\n", ADAPTOR, PASS_RATE);
	
	while (read_fastq(p_read) != NULL )
	{	
		num_read++;
		
		if ( adaptor_remove(p_read) )
		{
			contain_adaptors++;
            
            if ( strlen(p_read->sequence) >= min_length_pass_filter) print_fastq(p_read);
            else too_short++;
		}
	}
	
	
	fprintf(stderr, "Total number of reads:			 %ld\n", num_read);
	fprintf(stderr, "number of reads have adaptor :       %ld, that is %ld%% of total reads\n", contain_adaptors, contain_adaptors*100/num_read);
    fprintf(stderr, "number of reads have adaptor and longer than %d :       %ld, that is %ld%% of reads contail adaptors\n", min_length_pass_filter, contain_adaptors - too_short, (contain_adaptors - too_short)*100/contain_adaptors);
    
    fprintf(stdin, "done!\n");
    
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

void print_fastq(const struct fastq* read)
{
	printf("%s\n", read->name);
	printf("%s\n", read->sequence);
	printf("%s\n", read->tag);
	printf("%s\n", read->score);
}


Bool	found_adaptors(char* read)
{
	const char*	key = ADAPTOR;
    float		total=0, match=0;
    float   rate;

	while ( *read && *key)
	{
		total++;
        if ( *read == *key) match++;
        rate = match / total;
        if ( (total >= MIN_MATCH) && (rate >= PASS_RATE)) return TRUE;
        
        read++;
        key++;
	}
    return FALSE;
}
		
Bool	adaptor_remove( struct fastq* read)
{
	char*	p = read->sequence;
	char	*p_temp;
	

	while ( *p )  // go over the whole sequence string
	{
        if ( found_adaptors(p) )
        {
            *p = '\0';
            
            p_temp = p - read->sequence + read->score;
            *p_temp = '\0';
            
            return TRUE;
        }
        
        p++;
	}
    
	return FALSE;

}






