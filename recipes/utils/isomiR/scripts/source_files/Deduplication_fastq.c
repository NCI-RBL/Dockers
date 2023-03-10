
 
 
 
 
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

struct node_fastq
{
	struct fastq			read;
	int					abundance; 
	struct node_fastq * 	next;
};

char*	read_fastq(struct fastq*);

void		print_fastq(const struct fastq*);

void 	print_seq(const struct fastq* );

int	search(const struct fastq*, const char* );   /* first argument is the pointer to read, second is pointer to searching keyword, return is how many times the motif was found */ 

Bool	found_adaptors(const struct fastq*, const char*, const char*, long int distance[]);	/* return TURE if two motifs are found, otherwise return FALSE, distance between two motif is stored in distance array */

Bool	adaptor_remove( struct fastq*, const char*, const char*, const int);  /* remove all sequence after the last of motif1, then fill tag */

Bool	is_redundant( struct node_fastq *, struct fastq *);

struct node_fastq * add_node ( struct node_fastq *, struct fastq * );  /* return pointer to the last node */

int print_all_node(struct node_fastq *);	/* return number of node printed */

int main(int argc, char* argv[])		/* argv[1] and argv[2] will be the pointers to search motifs */
{
	if (argc != 1) 
	{
		fprintf(stderr, "\n deduplication_fastq.exe < input.fastq \n"); 
		return 0;
	}

	long int 	num_read = 0, num_redundant = 0;
	long int		len[LEN_SEQ+1] = {0};
	int			length, i;

	

	struct fastq	read_temp;
	struct fastq 	*p_read = &read_temp;		/* p_read point to the current read */

	struct node_fastq		node_first;
	struct node_fastq 	*p_node_head = &node_first;		/* p_node point to the first node */
	struct node_fastq		*p_node_last = p_node_head;

	if ((read_fastq(p_read) != NULL )) 
		{
			node_first.read		 = *p_read;
			node_first.abundance   = 1;
			node_first.next		 = NULL;
		
			num_read++;
		}
							

	while (read_fastq(p_read) != NULL )
	{	
		num_read++;
							
		if (is_redundant(p_node_head, p_read)) 	num_redundant++;
		else			p_node_last = add_node ( p_node_last, p_read);		 	
	}
	
	fprintf(stderr, "\n\n\n");
	fprintf(stderr, "Total number of reads:	%15ld\n", num_read);
	fprintf(stderr, "Total number of redundant reads:%6ld that is %ld%% of total reads\n", num_redundant, num_redundant*100/num_read);
	
	print_all_node(p_node_head);

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

Bool	is_redundant( struct node_fastq *p_node_head, struct fastq *p_read)
{
	struct node_fastq * p;
	
	for ( p = p_node_head; p != NULL ; p = p->next)
	{
		if ( (strcmp((p->read).tag, p_read->tag) == 0 ) && (strcmp ((p->read).sequence, p_read->sequence) == 0))
		{
			p->abundance++;
			return TRUE;
		}
	}
	
	return FALSE;
} 
		

struct node_fastq * add_node ( struct node_fastq * p_last, struct fastq * p_read)  /* return pointer to the new last node */
{
	struct node_fastq *p_new;
	
	if( (p_new = malloc(sizeof(struct node_fastq))) == NULL) 
		{
			fprintf(stderr, "Not enough memory!\n");
			exit (1);
		}
	
	p_new->read = * p_read;
	p_new->abundance = 1;
	p_new->next = NULL;

	p_last->next = p_new;

	return p_new;	
}

int print_all_node(struct node_fastq * p_node_head)	/* return number of node printed */
{
	long int count = 0;
	struct node_fastq * p;

	for ( p = p_node_head; p != NULL ; p = p->next)
	{
		print_fastq(&(p->read));
		count++;			
	}

	return count;
}






