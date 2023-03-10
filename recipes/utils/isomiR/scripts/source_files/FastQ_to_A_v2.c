/* this grogram converts fastq file to fasta file  */

// v2 Jan-20-2015
 

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define Bool		int
#define	TRUE		1
#define	FAIL        0
#define PASS        1

#define	LEN_NAME	75
#define	LEN_SEQ		150
#define	LEN_TAG		15

#define OFFSET      33

// date structure is here........

struct fastq
{
	char        name[LEN_NAME+1];
	char        sequence[LEN_SEQ+1];
	char		tag[LEN_TAG+1];
	char		score[LEN_SEQ+1];
};

// Function prototypes are here......

struct fastq*	read_fastq(struct fastq*);

void	print_fastq(const struct fastq*);

void	print_fasta( const struct fastq*);

Bool    check_Qscore(struct fastq* read);

// Global varables are here...

int 	num_read = 0;
float     num_pass = 0;

int     Qscore_fliter, Qscore_mask;

// main start from here....

int main(int argc, char* argv[])		/* argv[1] and argv[2] will be the min_scores for passing filter and being unmarked by N */
{
	if (argc != 3)
	{
		fprintf(stderr, "\n fastQ_to_A_v2.exe min_Qscore_pass_filter min_Qscore_unmarked < input.fastq > output.fa (Qscore_unmask should be higher than Qscore_pass_filter)\n\n");
		return 0;
	}
    
    Qscore_fliter = atoi(argv[1]) + OFFSET;
    Qscore_mask = atoi(argv[2]) + OFFSET;
    
    if (Qscore_mask < Qscore_fliter)
    {
        fprintf(stderr, "\n Qscore_unmask should be equal or higher than Qscore_pass_filter\n");
        return 0;
    }
    
    struct fastq	read_temp;
	struct fastq 	*p_read = &read_temp;		/* p_read point to the current read */
    
	while (read_fastq(p_read) != NULL )
	{
		num_read++;
        if(check_Qscore(p_read) == PASS) {num_pass++; print_fasta(p_read);}
	}
	
	
	fprintf(stderr, "\n\n\n");
	fprintf(stderr, "Total number of reads:	%15d\nNumber of reads passing filter: %15f, that is %.2f%% of total reads\n", num_read, num_pass, num_pass*100/num_read);

	return 0;	
}

Bool    check_Qscore(struct fastq* read)
{
    int index;
    int length;
    
    length = strlen(read->sequence);
    
    for(index = 0; index < length; index++)
    {
        if(read->score[index] >= Qscore_mask)  continue;         // good score, higher than mask
        if(read->score[index] < Qscore_fliter) return FAIL;      // low score, discard this read
        read->sequence[index] = 'N';
    }
    
    return PASS;
}



struct fastq* read_fastq(struct fastq* read)
{
	if (gets(read->name) == NULL) return NULL;
	gets(read->sequence);
	gets(read->tag);
	gets(read->score);
	
	return read;
}

void print_fastq(const struct fastq* read)
{
	printf("%s\n", read->name);
	printf("%s\n", read->sequence);
	printf("%s\n", read->tag);
	printf("%s\n", read->score);
}

void print_fasta(const struct fastq* read)
{
	putchar('>');
	printf("%s\n", read->name);
	printf("%s\n", read->sequence);
}
