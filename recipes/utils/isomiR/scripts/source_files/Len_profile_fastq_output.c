/* this grogram analyze the length profile of a fastq file */

// July-27-2018 Clean up the file to make it clear..
// Added region of analysis - define as CONSTANT
 
// July-30-2018
// add output - certain length as defined by user into txt files.
// the output will be all reads containing duplications. user can use excel later to de-duplicate those and get unique reads.


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define	LEN_NAME	175
#define	LEN_SEQ		175
#define	LEN_TAG		25

#define MIN_LEN		6
#define MAX_LEN		23

struct fastq 
{
	char 	name[LEN_NAME+1];
	char 	sequence[LEN_SEQ+1];
	char	tag[LEN_TAG+1];
	char	score[LEN_SEQ+1];
};

char*	read_fastq(struct fastq*);

int main(int argc, char* argv[])    // one argument, read length.
{
	if (argc != 2)
	{
		fprintf(stderr, "\n len_profile_fastq.exe filter_length < input.fastq > output.txt\n");
		return 0;
	}

	long int 	num_read = 0;
	long int	len[MAX_LEN+1] = {0};
	int		length, i;
    int     filter_length;
    
	struct fastq	read_temp;
	struct fastq 	*p_read = &read_temp;	/* p_read point to the current read */
    
    filter_length = atoi(argv[1]);
	
	while (read_fastq(p_read) != NULL )
	{
		length = strlen(p_read->sequence);
        if (length == filter_length) fprintf(stdout, "%s\n", p_read->sequence);
	}
}

char* read_fastq(struct fastq* read)
{
	if (gets(read->name) == NULL) return NULL;
	gets(read->sequence);
	gets(read->tag);
	gets(read->score);
	
	return read->name;
}

