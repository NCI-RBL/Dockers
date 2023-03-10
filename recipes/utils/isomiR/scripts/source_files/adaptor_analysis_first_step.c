/* this grogram search specific sequence motifs of adaptors in Read 1 of a fastq file*/
/* and then print out the statistics */
/* search_motif MOTIF1 MOTIF2 < file.fastq */
/* searching for the last of motif1 and the first of motif2 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define Bool	int
#define	TRUE	1
#define	FALSE	0

#define	LEN_NAME	75
#define	LEN_SEQ		75

#define MIN_DIS		14
#define MAX_DIS		20
	
struct fastq 
{
	char 	name[LEN_NAME+1];
	char 	sequence[LEN_SEQ+1];
	char	third[2];
	char	score[LEN_SEQ+1];
};

char*	read_fastq(struct fastq*);

void	print_fastq(const struct fastq*);

void 	print_seq(const struct fastq* );

int	search(const struct fastq*, const char* );   /* first argument is the pointer to read, second is pointer to searching keyword, return is how many times the motif was found */ 

Bool	found_adaptors(const struct fastq*, const char*, const char*, long int distance[]);	/* return TURE if two motifs are found, otherwise return FALSE, distance between two motif is stored in distance array */

int main(int argc, char* argv[])		/* argv[1] and argv[2] will be the pointers to search motifs */
{
	if (argc != 3) 
	{
		printf("\nsearch_motif MOTIF1 MOTIF2 < file.fasta\n"); 
		return 0;
	}

	long int 	num_read = 0, discarded = 0, contain_adaptors = 0;
	long int	distance[LEN_SEQ] = {0}, dis_short = 0, dis_long = 0;
	int	i;

	struct fastq	read_temp;
	struct fastq 	*p_read = &read_temp;	/* p_read point to the current read */

	while (read_fastq(p_read) != NULL )
	{	
		num_read++;
		
		if (found_adaptors(p_read, argv[1], argv[2], distance)) 	contain_adaptors++;
		else								discarded++;
	}
	
	putchar('\n');
	printf("Total number of reads:			 %15ld\n", num_read);
	printf("number of reads have two adaptors:       %15ld, that is %ld%% of total reads\n", contain_adaptors, contain_adaptors*100/num_read);
	printf("number of reads have been discarded:     %15ld, that is %ld%% of total reads\n", discarded, discarded*100/num_read);

	for (i = 1; i <= MIN_DIS; i++) dis_short += distance[i];
	for (i = MAX_DIS; i < LEN_SEQ; i++) dis_long += distance[i];
	
	putchar('\n');
	printf("Distance between adaptors       Number of reads     Percentage\n");
	printf("       <=%d             %15ld        %10ld%%\n", MIN_DIS, dis_short, dis_short*100/num_read);
	for (i = MIN_DIS + 1; i < MAX_DIS; i++)
		printf("         %d             %15ld        %10ld%%\n", i, distance[i], distance[i]*100/num_read);
	printf("       >=%d             %15ld        %10ld%%\n", MAX_DIS, dis_long, dis_long*100/num_read);        	

	return 0;	
}

char* read_fastq(struct fastq* read)
{
	if (gets(read->name) == NULL) return NULL;
	gets(read->sequence);
	gets(read->third);
	gets(read->score);
	
	return read->name;
}

void print_fastq(const struct fastq* read)
{
	printf("%s\n", read->name);
	printf("%s\n", read->sequence);
	printf("%s\n", read->third);
	printf("%s\n", read->score);
}

void print_seq(const struct fastq* read)
{
	printf("%s\n", read->sequence);
}

int  search(const struct fastq* read, const char* keyword )
{
	const char*	p = read->sequence;
	int	count = 0;

	while ( (p = strstr(p, keyword)) != NULL)
	{
		count++;
		p++;
	}

	return count;
}

Bool	found_adaptors(const struct fastq* read, const char* motif1, const char* motif2, long int* distance)
{
	const char*	p = read->sequence;
	char	*p_motif1_start, *p_motif1_end, *p_motif2_start, *p_temp;
	int		len_motif1, dis;
	
	len_motif1 = strlen(motif1);

	while ( (p_temp = strstr(p, motif1)) != NULL) 
	{
		p_motif1_start = p_temp;
		p = p_motif1_start+1; 
	}	

	if (p == read->sequence) return FALSE;

	p_motif1_end = p_motif1_start + len_motif1 ;		/* point to the next one after the end of motif1 */

	if (p_motif1_end >= read->sequence + LEN_SEQ) return FALSE;
	
	if ((p_motif2_start = strstr( p_motif1_end, motif2)) == NULL) return FALSE;
		 		
	dis = p_motif2_start - p_motif1_end;

	distance[dis]++;

	return TRUE;

}
		







