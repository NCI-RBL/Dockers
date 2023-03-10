/********************************************************************************************************************************************
 *                                                                                                                                          *
 * this grogram deal with raw data in the format of fastq, assuming the read is read1 which contains 15nt randam sequences tag at the end   *
 * of 3'end                                                                                                                                 *
 * program reads fasta file of read1 from stdin. The final output is sequnences ready to be mapped in fastq, which will be print to stdout  *
 * program will also generate information regarding the sequencing processing, print to stderr                                              *
 *                                                                                                                                          *
 *             Commond --    analyze_read1 mismatch_allowed minimal_size_pass_filter < input.fastq > output.fastq                           *
 *                                                                                                                                          *
 *              Constants defined are   1.  length of the read  -- 150                                                                      *
 *                                      2.  length of tag       -- 15                                                                       *
 *                                      3.  adaptor motifs info -- CTGAC-15N-TGGAATTCTC (30nt tag with 20nt recognable sequence)            *
 *                                      4.  allow mismatches- check pattern CTGAC-15N-TGGAATTCTC                                            *
 *                                      5.  loose de-dup condition : a.                                                      *
 ********************************************************************************************************************************************/

//  2014-Dec-18 : v2
//  find the first CTGAC-15N-TGGAA within mismatch allowrance
//  shortage for motif2 is allowed but counted in mismatch allowrance

//  2015-Jan-30 :   v3
//  changed the file name to trim_read1_miv_v3.c
//  this one changed the way we looking for adaptors -- instead of counting mismatch allowrance, using mininal matches and pass_rate
//  the inputs are changed accordingly

//  2015-Mar-13 :   v4
//  adaptor sequence changed to match the 3' adaptor used. 8nt (CTGTTAAC) + 15N
//  changed the way we are looking for mismatches. using pass-rate and minimal number of match.
//  Both 3' adaptor searching, 5' adaptor filtering and 15N plus read comparing not are now using pass_rate and minimal number of matches.
//  Using the first part of read sequences as key, instead of the 15N. More trustable at the 5' end? Tried.. can capture 2~3 fold more of the
//  duplicated reads. 3.7k -> 8.3k. However, consume too much time. So changed back!

//  2017-Aug-17:    v5
//  make the read length longer so it can analyze longer reads!

//  2017-Oct-16:    not passed
//  the output will be those not passed.

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define Bool		int
#define	TRUE		1
#define	FALSE       0
#define SUCCESS     1
#define FAIL        0

#define	LEN_NAME	150              /* need to define the data structure for fastq read */
#define	LEN_SEQ		150
#define	LEN_TAG		16

#define ADAPTOR		"CTGTTAACNNNNNNNNNNNNNNNTGGAATTCTCGGGT"

#define KAY_SIZE            7
#define HASH_TABLE_SIZE     16384        //key word size is 4^7


#define THREE_ADAPTOR_PASS_RATE         0.909
#define THREE_ADAPTOR_MIN_MATCH         8           // the CTGTTAAC are match OR one mismatches need to be compensated by three matches after 15N
#define FIVE_ADAPTOR_PASS_RATE          0.8
#define FIVE_ADAPTOR_MIN_MATCH          15
#define TAG_PASS_RATE                   0.8         // allowing up to 2 mimatches in the 15N tag..
#define TAG_MIN_MATCH                   13
#define READ_PASS_RATE                  0.8
#define READ_MIN_MATCH                  10

#define FIVE_ADAPTOR    "CTACACGACGCTCTTCCGATCT"    //5' adaptor seq



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

struct      node_fastq *   hash_table[ HASH_TABLE_SIZE ] = {NULL};       // hash table store the pointer to the first node contain the same key word


/****************************************/


/*                                         FUNTION PROTOTYPE HERE                          */

/* GENERAL IO FUNCTIONS */

struct node_fastq*       read_fastq(void);     /* read from stdin, created a node, store seq in, return the pointer to node if sucessful when reading 
                                                fails, return a NULL pointer */

void		print_fastq(struct fastq*);

int         print_all_nodes(void);	// print all nodes based from the hash table. return number of reads printed.



/* ADAPTOR SEARCH AND REMOVAL */

Bool        remove_adaptor( struct node_fastq * );    /* search then remove adaptor, turn status -- sucess or fail */

Bool        found_adaptors(char* read);

Bool        compare_str_rate(char* read, char* key, float pass_rate, int min_match);

/* DE-DUPLICATION */

Bool        is_duplicated( struct node_fastq*);

long int    key_value(char*);

Bool        is_node_same( struct node_fastq *p_node_head, struct node_fastq *p_read);





/************************************************************************************************************************************/


int main(int argc, char* argv[])		/* One argument allowed */
{
	if (argc != 2)
	{
		fprintf(stderr, "\n trim_read1_v4.exe min_length_pass_filter < input.fastq > output.fastq \n");
		return 0;
	}
	int         min_length_pass_filter = atoi (argv[1]);
	
    int			length;
	struct node_fastq* p_cur;
    int         num_total_reads = 0, num_contain_adaptor = 0, num_empty = 0, num_too_short = 0, num_5adapotr_contamination = 0, num_final = 0, num_duplicated = 0, num_no_adaptor=0;
    
	while ((p_cur = read_fastq()) != NULL )
	{	
		num_total_reads++;
        
        if (!remove_adaptor(p_cur)) {print_fastq (&(p_cur->read));  num_no_adaptor++;}
       
    }
    
    fprintf(stderr, "\nTotal number of reads:             %15d\n", num_total_reads);
    fprintf(stderr, "Number of reads have no adaptor :       %15d, that is %d%% of total reads\n", num_no_adaptor, num_no_adaptor*100/num_total_reads);
        
 /*       if (remove_adaptor(p_cur))
        {
            num_contain_adaptor++;
            
            length = strlen(p_cur->read.sequence);
            
            if (length == 0) num_empty++;
            else if (length <= min_length_pass_filter) num_too_short++;
            else if (compare_str_rate(p_cur->read.sequence, FIVE_ADAPTOR, FIVE_ADAPTOR_PASS_RATE, FIVE_ADAPTOR_MIN_MATCH)) num_5adapotr_contamination++;
            else if (is_duplicated(p_cur))  num_duplicated++;
        }
	}

	num_final = print_all_nodes();
    
    fprintf(stderr, "\nTotal number of reads:			 %15d\n", num_total_reads);
    fprintf(stderr, "Number of reads have adaptor :       %15d, that is %d%% of total reads\n", num_contain_adaptor, num_contain_adaptor*100/num_total_reads);
    fprintf(stderr, "Number of empty reads :       %15d, that is %d%% of reads passing adaptor removal\n", num_empty, num_empty*100/num_contain_adaptor);
    fprintf(stderr, "Number of reads shorter than %d :       %15d, that is %d%% of reads passing adaptor removal\n", min_length_pass_filter, num_too_short, num_too_short*100/num_contain_adaptor);
    fprintf(stderr, "Number of 5' adaptor contamination reads :       %15d, that is %d%% of reads passing adaptor removal\n", num_5adapotr_contamination, num_5adapotr_contamination*100/num_contain_adaptor);
    fprintf(stderr, "Number of duplicated reads :       %15d, that is %d%% of reads passing adaptor removal\n", num_duplicated, num_duplicated*100/num_contain_adaptor);
    fprintf(stderr, "Number of final reads :       %15d, that is %d%% of reads passing adaptor removal\n", num_final, num_final*100/num_contain_adaptor);
 */
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

void print_fastq(struct fastq* read)
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

int print_all_nodes(void)           // return number of reads printed...
{
    int    index;
    int    num_reads = 0;
    struct node_fastq * p_nod;
    struct fastq      * p_read;
    int    length;
    
    for ( index = 0; index < HASH_TABLE_SIZE; index++)
    {
        for ( p_nod = hash_table[index]; p_nod != NULL; p_nod = p_nod->next)
        {
            p_read = & (p_nod->read);
            print_fastq(p_read);
            num_reads++;
            
            
        }
    }
    return num_reads;
}

/* ADAPTOR SEARCH AND REMOVAL */


Bool	found_adaptors(char* read)
{
                char*       key = ADAPTOR;
                int         count = 0;
                float		total=0, match=0;
                float       rate;
                
                while ( *read && *key)
                {
                    total++;
                    if ( *read == *key) match++;
                    rate = match / total;
                    if ( (match >= THREE_ADAPTOR_MIN_MATCH) && (rate >= THREE_ADAPTOR_PASS_RATE)) return TRUE;
                    
                    if (++count == 8) {read+=15; key+=15;}
                    else {read++; key++;}
                }
                return FALSE;
            }
                            
Bool	remove_adaptor( struct node_fastq* p_nod)
{
        char	*p = p_nod->read.sequence;
        struct fastq    *read = &(p_nod->read);
        char	*p_temp;
                
                
        while ( *p )  // go over the whole sequence string
        {
            if ( found_adaptors(p) )
            {
                *p = '\0';                                      //mark the end of the sequence
                        
                p_temp = p - read->sequence + read->score;
                *p_temp = '\0';                                 // this two lines set the score string length the same as sequence
                        
                *(p + 23) = '\0';                               //mark the end of the 15N
                strcpy(p_nod->read.tag, p + 8);               // moved tag seq to third line
                
                return TRUE;
            }
                    
            p++;
        }
                
        return FALSE;
}


/* DE-DUPLICATION */

Bool    is_duplicated( struct node_fastq* p_nod)
{
    long int index;
    struct node_fastq * p_cur;
    
    
    index = key_value(p_nod->read.tag);
    
    if (hash_table[index] == NULL)
    {
        hash_table[index] = p_nod;
        return FALSE;
    }
    else
    {
        for (p_cur = hash_table[index]; p_cur != NULL; p_cur=p_cur->next)
        {
            if (is_node_same(p_cur, p_nod))  return TRUE;
        }
    }
    
    p_nod->next = hash_table[index];
    hash_table[index] = p_nod;
    return FALSE;
}



Bool    is_node_same( struct node_fastq *p_node_1, struct node_fastq *p_node_2)
{
	
    if ( compare_str_rate(p_node_1->read.tag, p_node_2->read.tag, TAG_PASS_RATE, TAG_MIN_MATCH) && compare_str_rate(p_node_1->read.sequence, p_node_2->read.sequence, READ_PASS_RATE, READ_MIN_MATCH)) return TRUE;
	else    return FALSE;
				
} 
		
long int    key_value(char* p_key)              /* caculate key value for a key of KAY_SIZE nt tag sequences */
{
    long int    value = 0;
    int         i;
    
    for (i=0; i < KAY_SIZE ; i++)
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
                            
Bool	compare_str_rate(char* read, char* key, float pass_rate, int min_match)
{
    float		total=0, match=0;
    float       rate;
    
    while ( *read && *key)
    {
        total++;
        if ( *read == *key) match++;
        rate = match / total;
        if ( (match >= min_match) && (rate >= pass_rate)) return TRUE;
                    
        read++;
        key++;
    }
    return FALSE;
}

