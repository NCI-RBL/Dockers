/* this grogram removes adaptors in Read 1 of a fastq file */
/* the tag sequence (15nt) will be stored in third line of each read */
/* and then print out the statistics */
/* adaptor_remove MOTIF1 MOTIF2 < file.fastq > output.fastq */
/* searching for the last of motif1 and the first of motif2 */
/* for this version: motif1 is CTGAC and motif2 is GAATTCT  */ 

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define Bool	int
#define	TRUE	1
#define	FALSE	0

#define	LEN_NAME	75
#define	LEN_SEQ		75
#define	LEN_TAG		15

#define MIN_DIS		14
#define MAX_DIS		20
#define DIS_MOTIFS	17

#define MOTIF1		"CTGAC"
#define	MOTIF2		"GAATTCT"
	
struct fastq 
{
	char 	name[LEN_NAME+1];
	char 	sequence[LEN_SEQ+1];
	char	tag[LEN_TAG+1];
	char	score[LEN_SEQ+1];
};

char*	read_fastq(struct fastq*);

void	print_fastq(const struct fastq*);

void 	print_seq(const struct fastq* );

int	search(const struct fastq*, const char* );   /* first argument is the pointer to read, second is pointer to searching keyword, return is how many times the motif was found */ 

Bool	found_adaptors(const struct fastq*, const char*, const char*, long int distance[]);	/* return TURE if two motifs are found, otherwise return FALSE, distance between two motif is stored in distance array */

Bool	adaptor_remove( struct fastq*, const char*, const char*, const int);  /* remove all sequence after the last of motif1, then fill tag */

int main(int argc, char* argv[])		/* argv[1] and argv[2] will be the pointers to search motifs */
{
	if (argc != 1) 
	{
		printf("\nadaptor_remove  < input.fastq > output.fastq\n"); 
		return 0;
	}

	long int 	num_read = 0, discarded = 0, contain_adaptors = 0;
	long int	distance[LEN_SEQ] = {0}, dis_short = 0, dis_long = 0;
	int	i;

	struct fastq	read_temp;
	struct fastq 	*p_read = &read_temp;	/* p_read point to the current read */

	fprintf(stderr, "\nMotif one is %s and Motif two is %s\n", MOTIF1, MOTIF2);
	
	while (read_fastq(p_read) != NULL )
	{	
		num_read++;
		
		if ( adaptor_remove(p_read, MOTIF1, MOTIF2, DIS_MOTIFS) ) 	
		{
			contain_adaptors++;
			print_fastq(p_read);
		}
		else	discarded++;
	}
	
	
	fprintf(stderr, "Total number of reads:			 %15ld\n", num_read);
	fprintf(stderr, "number of reads have two adaptor :       %15ld, that is %ld%% of total reads\n", contain_adaptors, contain_adaptors*100/num_read);
	fprintf(stderr, "number of reads have been discarded:     %15ld, that is %ld%% of total reads\n", discarded, discarded*100/num_read);

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
		
Bool	adaptor_remove( struct fastq* read, const char* motif1, const char* motif2, const int dis_motifs)
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
		 		
	if (dis_motifs != (p_motif2_start - p_motif1_end)) return FALSE;

	*(p_motif2_start - 2) = '\0';     

	strcpy(read->tag, p_motif1_end);

	*p_motif1_end = '\0';

	return TRUE;

}






