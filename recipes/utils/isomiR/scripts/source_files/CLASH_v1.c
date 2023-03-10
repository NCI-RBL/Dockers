
// 2015 Jan 6th  Rewrite the old version to analyze data from illumina small RNA kit.

// search 3' adaptor sequnce and remove it.

// based on the percentage of match- for example 0.8 means 1 mismatch is allowed out of 5nt

// find the first one then remove the rest.

// maybe I need to add a feature to remove all contaminat

// maybe I need to give stat about length distribution

// 2016 Jan Find out that the first 5nt of hsa-miR-1 is excact the same as adaptor. Changed the pass-rate to 0.9 and minimal length to 10. 

// 2016 Jan 20th, change the minimal length of reads to 10nt -> to analyze the trimmed reads of miR-27a/b/s

// 2017 Oct 20th, change the length of reads to 110

// 2017 Oct 20th, change the file to analyze the clash results from Darnell group..

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define Bool	int
#define	TRUE	1
#define	FALSE	0
#define SUCCESS     1
#define FAIL        0

#define	LEN_NAME	110
#define	LEN_SEQ		110
#define	LEN_TAG		110

#define ADAPTOR		"TGGAATTCTCTCGGGTGCCAAGG"
#define	PASS_RATE           0.9
#define MIN_MATCH           10

#define MIR9C      "TCTTTGGTTATC"
#define MIR9A      "CTTTGGTTATCT"

#define MIN_LEN             10
	
struct fastq 
{
	char 	name[LEN_NAME+1];
	char 	sequence[LEN_SEQ+1];
	char	tag[LEN_TAG+1];
	char	score[LEN_SEQ+1];
};

Bool    read_fastq(struct fastq *, FILE*);

void	print_fastq(const struct fastq*);

Bool	found_adaptors(char*);	/* return TURE if adaptor is found, otherwise return FALSE */

Bool	adaptor_remove( struct fastq* );  // removed return TRUE, otherwise return FALSE

int main(int argc, char* argv[])		
{
	if (argc != 3)
	{
		printf("\n CLASHV_v1.exe  input.fastq  output\n");
		return 0;
	}

    FILE                *fp_input, *fp_output_miR9C, *fp_output_miR9A;
    char                miR9C[200], miR9A[200];
    
    
    if ( (fp_input = fopen(argv[1], "r")) == NULL )                     // open motif file, making sure it is sucessful..
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[1]);
        exit (1);
    }
    
    strcpy (miR9C, argv[2]);
    strcat (miR9C, "_miR9-C.fastq");
    strcpy (miR9A, argv[2]);
    strcat (miR9A, "_miR9-A.fastq");
    
    
    fp_output_miR9C = fopen( miR9C, "w");
    fp_output_miR9A = fopen( miR9A, "w");
    
	long int 	num_read = 0,  contain_miR9C = 0, contain_miR9A = 0;
    int         size;

	struct fastq	read_temp;
	struct fastq 	*p_read = &read_temp;	/* p_read point to the current read */
    
    char    miRNA[30];
    char    first_section[60];
    char    target_section[60];

	fprintf(stderr, "\n miR9-Canonical sequence is %s and miR9-Canonical sequence is %s\n", MIR9C,MIR9A);
	
	while (read_fastq(p_read, fp_input) != FAIL )
	{	
		num_read++;
        strncpy(miRNA, p_read->sequence+25, 26);
        strncpy(first_section, p_read->sequence+20, 31);
        strcpy(target_section, p_read->sequence+51);
        
		if (strncmp(miRNA, MIR9C, 12) == 0)
        {
            size = strlen (target_section);
            p_read->score[size] = '\0';
            
            fprintf(fp_output_miR9C, "%s", p_read->name);
            fprintf(fp_output_miR9C, "%s", target_section);
            fprintf(fp_output_miR9C, "+ %s\n", first_section);
            fprintf(fp_output_miR9C, "%s\n", p_read->score);
            
            contain_miR9C++;
            continue;
        }
        
        if (strncmp(miRNA, MIR9A, 12) == 0)
        {
            size = strlen (target_section);
            p_read->score[size] = '\0';
            
            fprintf(fp_output_miR9A, "%s", p_read->name);
            fprintf(fp_output_miR9A, "%s", target_section);
            fprintf(fp_output_miR9A, "+ %s\n", first_section);
            fprintf(fp_output_miR9A, "%s\n", p_read->score);
            
            contain_miR9A++;
            continue;
        }
	}
	
	
	fprintf(stderr, "Total number of reads:			 %15ld\n", num_read);
	fprintf(stderr, "number of reads contain miR9-C :       %15ld, that is %ld%% of total reads\n", contain_miR9C, contain_miR9C*100/num_read);
    fprintf(stderr, "number of reads contain miR9-A :       %15ld, that is %ld%% of total reads\n", contain_miR9A, contain_miR9A*100/num_read);
    
	return 0;	
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

Bool    read_fastq(struct fastq* p_read, FILE* fp_input)
{
    if (fgets(p_read->name, LEN_SEQ, fp_input) == NULL) return FAIL;
    if (fgets(p_read->sequence, LEN_SEQ, fp_input) == NULL) return FAIL;
    if (fgets(p_read->tag, LEN_SEQ, fp_input) == NULL) return FAIL;
    if (fgets(p_read->score, LEN_SEQ, fp_input) == NULL) return FAIL;
    
    return SUCCESS;
}





