/* this grogram converts fastq file to fasta file  */
 

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define 	Bool		int
#define	TRUE		1
#define	FALSE	0

#define	LEN_NAME		75
#define	LEN_SEQ		75
#define	LEN_TAG		15

#define MIN_DIS		14
#define MAX_DIS		20
#define DIS_MOTIFS	17

#define MOTIF1		"CTGAC"
#define MOTIF2		"GAATTCT"
	
struct fastq 
{
	char 	name[LEN_NAME+1];
	char 	sequence[LEN_SEQ+1];
	char		tag[LEN_TAG+1];
	char		score[LEN_SEQ+1];
};

char*	read_fastq(struct fastq*);

void	print_fastq(const struct fastq*);

void 	print_seq(const struct fastq* );

void	print_fasta( const struct fastq*);


int main(int argc, char* argv[])		/* argv[1] and argv[2] will be the pointers to search motifs */
{
	if (argc != 1) 
	{
		fprintf(stderr, "\n fastq-fasta.exe < input.fastq > output.fa\n"); 
		return 0;
	}

	long int 	num_read = 0;

	struct fastq	read_temp;
	struct fastq 	*p_read = &read_temp;		/* p_read point to the current read */

	while (read_fastq(p_read) != NULL )
	{	
		num_read++;
							
		print_fasta(p_read);	 	
	}
	
	fprintf(stderr, "\n\n\n");
	fprintf(stderr, "Total number of reads:	%15ld\n", num_read);

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

void print_seq(const struct fastq* read)
{
	printf("%s\n", read->sequence);
}

void print_fasta(const struct fastq* read)
{
	putchar('>');
	printf("%s\n", read->name);
	printf("%s\n", read->sequence);
}
