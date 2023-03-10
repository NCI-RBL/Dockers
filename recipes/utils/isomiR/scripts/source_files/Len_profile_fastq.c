/* this grogram analyze the length profile of a fastq file */

// July-27-2018 Clean up the file to make it clear..
// Added region of analysis - define as CONSTANT
 

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define	LEN_NAME	175
#define	LEN_SEQ		175
#define	LEN_TAG		25

#define MIN_LEN		15
#define MAX_LEN		70

struct fastq 
{
	char 	name[LEN_NAME+1];
	char 	sequence[LEN_SEQ+1];
	char	tag[LEN_TAG+1];
	char	score[LEN_SEQ+1];
};

char*	read_fastq(struct fastq*);

int main(int argc, char* argv[])    // NO Argument
{
	if (argc != 1) 
	{
		fprintf(stderr, "\n len_profile_fastq.exe < input.fastq \n"); 
		return 0;
	}

	long int 	num_read = 0;
    long int    num_skipped = 0;
	long int	len[MAX_LEN+1] = {0};
	int		length, i;

	struct fastq	read_temp;
	struct fastq 	*p_read = &read_temp;	/* p_read point to the current read */
	
	while (read_fastq(p_read) != NULL )
	{	
		num_read++;
		
		length = strlen(p_read->sequence);
		if (length >= MIN_LEN && length <= MAX_LEN) len[length]++;
        else num_skipped++;
	}
	
	printf("\n");
    printf("Total number of reads:  \t%15ld\n", num_read);
    printf("Number of skipped reads:\t%15ld\n", num_skipped);
    printf("\n");
	printf("Size\tCounts\tPercentage\n");
	
	for (i = MIN_LEN; i <= MAX_LEN; i++)
		printf("%d\t%ld\t%.4f%%\n", i, len[i], (float)len[i]*100/(float)num_read);
    printf("\n");
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

