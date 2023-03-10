/********************************************************************************************************************************************
 *                                                                                                                                          *
 * this grogram deal with raw data in the format of fastq, assuming the read is read1 which contains 15nt randam sequences tag at the end   *
 * of 3'end                                                                                                                                 *
 * program reads fasta file of read1 from stdin. The final output is sequnences ready to be mapped in fastq, which will be print to stdout  *
 * program will also generate information regarding the sequencing processing, print to stderr                                              *
 *                                                                                                                                          *
 *             Commond --    analyze_read1 < input.fastq > output.fastq                                                                     *
 *                                                                                                                                          *
 *              Constants defined are   1.  length of the read  -- 75                                                                       *
 *                                      2.  length of tag       -- 15                                                                       *
 *                                      3.  adaptor motifs info -- sequences and the distance between these two                             *
 *                                      4.  Minimal Size passing filter -- 18                                                               *
 *                                                                                                                                          *
 ********************************************************************************************************************************************/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define Bool		int
#define	TRUE		1
#define	FALSE       0
#define SUCCESS      1
#define FAIL        0

#define	LEN_NAME	76              /* need to define the data structure for fastq read */
#define	LEN_SEQ		100
#define	LEN_TAG		16

#define MOTIF1		"CTGAC"         /* needed to identify the 3' adaptor sequence */
#define MOTIF2		"GAATTCT"
#define DIS_MOTIFS	17

#define HASH_TABLE_SIZE     1048576     /* key word size is TEN. 4^10 = 1048576 */

#define NUM_DIS     5               /* used in analysis -- display how many top ones */

#define MIN_READ_SIZE   18          /* any size smaller than this will not pass filter */

/*   DATA STRUCTRUE DEFINITITION HERE   */
	
struct fastq 
{
	char        name[LEN_NAME+1];
	char        sequence[LEN_SEQ+1];
	char		tag[LEN_TAG+1];
	char		score[LEN_SEQ+1];
};

struct node_fastq
{
	struct fastq			read;
	int                     abundance;
	struct node_fastq * 	next;
};

/****************************************/

/*   GOLOBAL VIRABLES DEFINED HERE   */

struct      node_fastq *   hash_table[ HASH_TABLE_SIZE ] = {NULL};       /* hash table store the pointer to the first node contain the same key word */

long int    num_total_reads = 0;
long int 	num_reads_missing_motif1 = 0, num_reads_missing_motif2 = 0, num_reads_wrong_motif_distance = 0;
long int    num_reads_pass_adaptor_remove = 0, num_reads_final = 0;

long int	len_before_deduplication[LEN_SEQ+1] = {0};                  /* used for length profile analysis */
long int	len_after_deduplication[LEN_SEQ+1]  = {0};
int         len_max_before[NUM_DIS+1] = {0};                            /* the max table store INDEX of len table */
int         len_max_after[NUM_DIS +1] = {0};

long int    num_reads_redundant = 0;                                   /* used for dedup analysis */
/****************************************/


/*                                         FUNTION PROTOTYPE HERE                          */

/* GENERAL IO FUNCTIONS */

struct node_fastq*       read_fastq(void);     /* read from stdin, created a node, store seq in, return the pointer to node if sucessful when reading 
                                                fails, return a NULL pointer */

void		print_fastq(const struct fastq*);

void        print_seq(const struct fastq* );

void        print_all_nodes(void);	/* print all nodes based from the hash table */



/* ADAPTOR SEARCH AND REMOVAL */

Bool        remove_adaptor( struct node_fastq * );  /* search then remove adaptor, turn status -- sucess or fail */



/* DE-DUPLICATION */

void        de_duplication( struct node_fastq*);

long int    key_value(char*);

Bool        is_node_same( struct node_fastq *p_node_head, struct node_fastq *p_read);


/* ANALYSIS */

void    find_max_lentable(long int *, int *);

void    print_sta_remove_adaptor(void);
void    print_sta_deduplication(void);



/************************************************************************************************************************************/


int main(int argc, char* argv[])		/* no arguments allowed */
{
	if (argc != 1) 
	{
		fprintf(stderr, "\n analyze_read1 < input.fastq > output.fastq \n");
		return 0;
	}

	
	int			length, i = 10;

	struct node_fastq* p_cur;

	while ((p_cur = read_fastq()) != NULL )
	{	
		num_total_reads++;
        
        if ( remove_adaptor(p_cur) == FAIL ) continue ;

        length = strlen( p_cur->read.sequence);
        len_before_deduplication[length]++;
        
        de_duplication(p_cur);   /* link p_cur into the hash table */
	
	}

	print_all_nodes();          /* this will also update the len_after_dup table while go throught the hash-table */

    find_max_lentable(len_before_deduplication, len_max_before);
    find_max_lentable(len_after_deduplication, len_max_after);
    

    print_sta_remove_adaptor();
    print_sta_deduplication();

	return 0;	
}

/***********************************************************************************************/

/* GENERAL IO FUNCTIONS */

struct node_fastq*       read_fastq(void)     /* read from stdin, created a node, store seq in, return the pointer to node if sucessful, when reading
                                                fails, return a NULL pointer */
{
	struct  node_fastq  *p_new;
    char                *pos;
    
	if( (p_new = malloc(sizeof(struct node_fastq))) == NULL)
    {
        fprintf(stderr, "Not enough memory!\n");
        exit (1);
    }
    
    if (fgets(p_new->read.name, LEN_NAME + 1, stdin) == NULL) return NULL;
    if ((pos=strchr(p_new->read.name, '\n')) != NULL)           *pos = '\0';
	
    fgets(p_new->read.sequence, LEN_SEQ +1, stdin);
	if ((pos=strchr(p_new->read.sequence, '\n')) != NULL)       *pos = '\0';
    
    fgets(p_new->read.tag, LEN_TAG +1, stdin);
	if ((pos=strchr(p_new->read.tag, '\n')) != NULL)            *pos = '\0';

    fgets(p_new->read.score, LEN_SEQ +1, stdin);
	if ((pos=strchr(p_new->read.score, '\n')) != NULL)          *pos = '\0';
	
	p_new->abundance = 1;
	p_new->next = NULL;
    
	return p_new;
}

void print_fastq(const struct fastq* read)
{
	int length;
    char* pos;
    
    length = strlen(read->sequence);
    
    printf("%s\n", read->name);
	printf("%s\n", read->sequence);
	
    printf("%c ", '+');
    printf("%s\n", read->tag);
	
    pos = read->score + length;
    *pos = '\0';
    printf("%s\n", read->score);
}

void print_seq(const struct fastq* read)
{
	printf("%s\n", read->sequence);
}


void print_all_nodes(void)
{
    long int    index;
    struct node_fastq * p_nod;
    struct fastq      * p_read;
    int         length;
    
    for ( index = 0; index < HASH_TABLE_SIZE; index++)
    {
        for ( p_nod = hash_table[index]; p_nod != NULL; p_nod = p_nod->next)
        {
            p_read = & (p_nod->read);
            length = strlen( p_nod->read.sequence);
            len_after_deduplication[length]++;
            if ( length >= MIN_READ_SIZE )
            {
                print_fastq(p_read);
                num_reads_final++;
            }
        }
    }
}

/* ADAPTOR SEARCH AND REMOVAL */

Bool        remove_adaptor( struct node_fastq * p_nod)  /* search then remove adaptor, turn status -- sucess or fail */
{
    char    *p_seq = p_nod->read.sequence;
	char	*p_motif1_start, *p_motif1_end, *p_motif2_start, *p_temp;
	int		len_motif1, dis;
	
	len_motif1 = strlen(MOTIF1);
    
	while ( (p_temp = strstr(p_seq, MOTIF1)) != NULL)   /* this loop find the LAST motif1 in read */
	{
		p_motif1_start = p_temp;
		p_seq = p_motif1_start+1;
	}
    
	if (p_seq == p_nod->read.sequence) { num_reads_missing_motif1++; return FAIL;}  /* cannot find motif1 */
    
	p_motif1_end = p_motif1_start + len_motif1 ;		/* point to the next one after the end of motif1 */
    
    if (p_motif1_end >= p_nod->read.sequence + LEN_SEQ) {num_reads_missing_motif2++; return FAIL;}
     
    if ((p_motif2_start = strstr( p_motif1_end, MOTIF2)) == NULL) {num_reads_missing_motif2++; return FAIL;} /* cannot find motif2 */
     
    if (DIS_MOTIFS != (p_motif2_start - p_motif1_end)) {num_reads_wrong_motif_distance++; return FAIL;} /* wrong spacing between motifs */
     
    *(p_motif2_start - 2) = '\0';
     
    strcpy(p_nod->read.tag, p_motif1_end);  /* moved tag seq to third line */
     
    *p_motif1_start = '\0';                 /* adaptor finally removed */
     
    num_reads_pass_adaptor_remove++;
    
    return SUCCESS;
}



/* DE-DUPLICATION */

void    de_duplication( struct node_fastq* p_nod)
{
    long int index;
    struct node_fastq * p_cur;
    
    
    index = key_value(p_nod->read.tag);
    
    if (hash_table[index] == NULL) hash_table[index] = p_nod;
    else
    {
        for (p_cur = hash_table[index]; p_cur != NULL; p_cur=p_cur->next)
        {
            if (is_node_same(p_cur, p_nod)) { p_cur->abundance++; num_reads_redundant++; break;}
        }
        
        if (p_cur == NULL)
        {
            p_nod->next = hash_table[index];
            hash_table[index] = p_nod;
        }
    }
        
    
}



Bool    is_node_same( struct node_fastq *p_node_1, struct node_fastq *p_node_2)         /* currently without allowing any mismatches */
{
	
    if ( (strcmp(p_node_1->read.tag, p_node_2->read.tag) == 0 ) && (strcmp(p_node_1->read.sequence, p_node_2->read.sequence) == 0)) return TRUE;
	else    return FALSE;
				
} 
		
long int    key_value(char* p_key)              /* caculate key value for a key of 10 nt tag sequences */
{
    long int    value = 0;
    int         i;
    
    for (i=0; i < 10; i++)
    {
        value *= 4;
        switch ( p_key[i] )
        {
            case 'A':   value+= 0; break;
            case 'T':   value+= 1; break;
            case 'G':   value+= 2; break;
            case 'C':   value+= 3; break;
        }
    }
    return  value;
}

/* ANALYSIS */

void    find_max_lentable(long int * table, int * max_table)           /* max_table store the INDEX of len table */
{
    int                 i, j, k, max_index;
    Bool                exist = FALSE;
    
    
    for (i = 1; i <= NUM_DIS; i++)
    {
        max_index = 1;
        
        for ( j =0; j <= LEN_SEQ; j++, exist = FALSE)
        {
            for (k = 1; k < i ; k++)
            {
                if (j == max_table[k]) exist = TRUE ;
            }
            
            if ( table[j] > table[max_index] && exist == FALSE ) max_index = j;
        }
        
        max_table[i] = max_index;
    }
}

void    print_sta_remove_adaptor(void)
{
    int i;
    
    fprintf(stderr, "\nTotal number of reads is :                                                   %ld\n", num_total_reads);
    fprintf(stderr, "\nNumber of reads missing the frist motif -%s is :                        %ld, that is %ld%% of total reads\n", MOTIF1, num_reads_missing_motif1, num_reads_missing_motif1*100/num_total_reads);
    fprintf(stderr, "Number of reads have frist motif but missing the second motif -%s is:  %ld, that is %ld%% of total reads\n", MOTIF2, num_reads_missing_motif2, num_reads_missing_motif2*100/num_total_reads);
    fprintf(stderr, "Number of reads have both motifs but distance in between is off:               %ld, that is %ld%% of total reads\n", num_reads_wrong_motif_distance , num_reads_wrong_motif_distance*100/num_total_reads);
    fprintf(stderr, "\nOverall, Number of reads passing adaptor-removal is:                         %ld, that is %ld%% of total reads\n", num_reads_pass_adaptor_remove , num_reads_pass_adaptor_remove*100/num_total_reads);
    
    
    fprintf(stderr, "\n\nLength profile after adaptor removal:\n\n");
    fprintf(stderr, "Most dominant lengthes are:\n");
    for (i = 1; i <= NUM_DIS; i++)
    {
        fprintf(stderr, "Top %d:\t\t %2d nt\t\t%10ld reads\t\t%2ld%% of total reads\n", i, len_max_before[i], len_before_deduplication[len_max_before[i]], len_before_deduplication[len_max_before[i]]*100/num_total_reads );
    }

    fprintf(stderr, "\n\nDetailed distribution:\n\n");
    
    fprintf(stderr, "Size\tNumber of reads\tPercentage\n");
	
	for (i = 0; i <= LEN_SEQ; i++)
		fprintf(stderr, "%d\t%8ld\t%ld%%\n", i, len_before_deduplication[i], len_before_deduplication[i]*100/num_total_reads);

}

void    print_sta_deduplication(void)
{
    int         i;
    long int index;
    long int count_null = 0;
    
    for (index = 0; index < HASH_TABLE_SIZE; index++ )
    {
        if ( hash_table[index] == NULL) count_null++;
    }
    
    fprintf(stderr, "\n\nNumber of non-redundant reads is:                         %ld, that is %ld%% of total reads\n", num_reads_pass_adaptor_remove - num_reads_redundant , 100 - num_reads_redundant*100/num_reads_pass_adaptor_remove);
    
    
    fprintf(stderr, "\n\n%ld%% of hash table is used\n\n", 100 - count_null*100/HASH_TABLE_SIZE);
    
    fprintf(stderr, "\n\nLength profile after de-deplication:\n\n");
    fprintf(stderr, "Most dominant lengthes are:\n");
    for (i = 1; i <= NUM_DIS; i++)
    {
        fprintf(stderr, "Top %d:\t\t %2d nt\t\t%10ld reads\t\t%2ld%% of total reads\n", i, len_max_after[i], len_after_deduplication[len_max_after[i]], len_after_deduplication[len_max_after[i]]*100/num_total_reads );
    }
    
    fprintf(stderr, "\n\nDetailed distribution:\n\n");
    
    fprintf(stderr, "Size\tNumber of reads\tPercentage\n");
	
	for (i = 0; i <= LEN_SEQ; i++)
		fprintf(stderr, "%d\t%8ld\t%ld%%\n", i, len_before_deduplication[i], len_before_deduplication[i]*100/num_total_reads);
    
    fprintf(stderr, "\n\nFinal reads passing size filter: %ld. That is %ld%% of reads passin adaptor removal and %ld%% of total reads\n\n", num_reads_final, num_reads_final*100/num_reads_pass_adaptor_remove, num_reads_final*100/num_total_reads);
    
}







































