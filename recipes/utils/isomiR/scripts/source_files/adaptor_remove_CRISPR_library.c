
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

void	print_fastq(const struct fastq*, FILE*);

Bool	found_adaptors(char*, char*);	/* return TURE if adaptor is found, otherwise return FALSE */

Bool	adaptor_remove( struct fastq* );  // removed return TRUE, otherwise return FALSE

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
    
    
	long int 	num_read = 0,  contain_adaptors = 0;

	struct fastq	read_temp;
	struct fastq 	*p_read = &read_temp;	/* p_read point to the current read */

	fprintf(stderr, "\n U6 Adaptor sequence is %s, tracrRNA Adaptor sequence is %s and passing_rate is %f\n", ADAPTOR_U6, ADAPTOR_tracrRNA, PASS_RATE);
	
	while (read_fastq(p_read) != NULL )
	{	
		num_read++;
		
		if ( adaptor_remove(p_read) )
		{
			contain_adaptors++;
            
            print_fastq(p_read, fp_output);
		}
        else print_fastq(p_read, fp_no_adaptor);
    }
	
	
	fprintf(stderr, "Total number of reads:			 %15ld\n", num_read);
	fprintf(stderr, "number of reads have adaptor :       %15ld, that is %ld%% of total reads\n", contain_adaptors, contain_adaptors*100/num_read);
    
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


Bool	found_adaptors(char* read, char* key)
{
	
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
	char    temp[21];           // To hold the sgRNA sequence.

	while ( *p )  // go over the whole sequence string
	{
        if ( found_adaptors(p, ADAPTOR_U6) && found_adaptors(p+40, ADAPTOR_tracrRNA))
        {
            
            p = p + 20;     // now point to the first nt of sgRNA
            *(p+20) = '\0'; // add null to the end of sgRNA
            
            strcpy (temp, p);
            strcpy (read->sequence, temp);
            
            *(read->score + 20) = '\0';   // trim the score to 20nt long
            
            return TRUE;
        }
        
        p++;
	}
    
	return FALSE;

}






