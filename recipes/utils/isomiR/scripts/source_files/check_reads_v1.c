
//  2015-Jan-7 : v1
//  This program take trimmed reads (xxx_ready.fastq) as input. Sorting out reads contain certain motifs and some statistics..
//  motifs are listed in a seperate file in fasta format. the name of this file is provided in the command line as the first argument.
//  check_reads_v1.exe motif_list input.fastq output.txt


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define Bool		int
#define	TRUE		1
#define	FALSE       0
#define SUCCESS     1
#define FAIL        0


#define	LEN_SEQ		100


#define LEN_MOTIF_NAME  50          // need to define the data structure of motif tree
#define LEN_MOTIF_SEQ   50

#define MIS_ALLOW       0           // used in searching motif in reads


/*   DATA STRUCTRUE DEFINITITION HERE   */

struct node_read                                        // used to hold each UNIQUE read
{
	char                    sequence[LEN_SEQ+1];
	int                     length;
    int                     abundance;
	struct node_read *      next;
};

struct node_motif                                       // used to hold each motif (sequence keyword)
{
    char                    name[LEN_MOTIF_NAME];
    char                    sequence[LEN_MOTIF_SEQ];
    int                     abundance;
    struct  node_read *     p_node_read;
    struct  node_motif *    next;
};


/****************************************/


/*                                         FUNTION PROTOTYPE HERE                          */

char*   strstr_mis (char* s1, char* s2, int n );        /* Returns a pointer to the first occurrence of str2 in str1, or a null pointer if str2 is not part of str1. */
                                                        /* n is the mismatch allowrance */
                                                        /* strstr_mis search the whole s1, and report the FIRST matching of s2 with the LOWEST mismatch number */
int     strcmp_mis ( char* s1, char* s2, int n);        /* s1 is the pointer to the subject string, s2 substring, n is the numbr of mismatches allowed */
                                                        /* for return, -1 - no match, mismatches are more than n, -2 - s1 is shorter then s2, m>=0 - remaining mismatches allowrance */
                                                        /* s2 need a \0 at the the end, no such requirment on s1 */


// ============================================================================================

struct node_motif*  read_motif(FILE* fp_motif);             //read one motif (two lines) from the file fp_motif pointing to. return a pointer if sucessful, otherwise return a NULL pointer.

char*               read_seq(char* read_cur, FILE* fp_input);               //read one sequence from fastq file.. return NULL when reaching end...

Bool        is_match(char* read, char* motif);             //currently no mismatch allowed.. return TRUE if match, otherwise return FAIL

struct node_read*        link_read(char* read, struct node_read* p_read_begin);  //link the read to the node_read tree


/************************************************************************************************************************************/


int main(int argc, char* argv[])		/* Three arguments allowed */
{
	
    FILE               *fp_motif, *fp_input, *fp_output;
    int                 num_total_reads = 0;
    struct node_motif*  p_motif_begin = NULL, *p_motif_cur;
    struct node_read*   p_read_cur;
    
    if (argc != 4)
	{
		fprintf(stderr, "\n check_reads_v1.exe motif_list input.fastq output.txt \n");
		return 0;
	}
    
    if ( (fp_motif = fopen(argv[1], "r")) == NULL )                     // open motif file, making sure it is sucessful..
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[1]);
        exit (1);
    }

	while ((p_motif_cur = read_motif(fp_motif)) != NULL )                       // start reading motif file to a node-tree
	{	
		p_motif_cur->next = p_motif_begin;
        p_motif_begin = p_motif_cur;
	}

    fclose(fp_motif);
    
    if ( p_motif_begin == NULL )                     // check if there are some motifs
    {
        fprintf(stderr, "\n %s file is empty!\n", argv[1]);
        exit (1);
    }

    if ( (fp_input = fopen(argv[2], "r")) == NULL )                     // open input file, making sure it is sucessful..
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[2]);
        exit (1);
    }

    char read_cur[LEN_SEQ + 1];
    
    while ( read_seq(read_cur, fp_input) != NULL)
    {
        num_total_reads++;
        
        for( p_motif_cur = p_motif_begin; p_motif_cur != NULL; p_motif_cur = p_motif_cur->next)
        {
            if ( is_match( read_cur, p_motif_cur->sequence))
            {
                p_motif_cur->abundance++;
                p_motif_cur->p_node_read = link_read(read_cur, p_motif_cur->p_node_read);
                break;
            }
        }
    }
    
    fclose(fp_input);
    
    fp_output = fopen(argv[3], "w");
    
    fprintf(fp_output, "Motif list file is %s; Input file is %s; Number of total reads is %d\n", argv[1], argv[2], num_total_reads);
    fprintf(fp_output, "Allowing no mismatches in motif searching!\n\n\n");
	
    for( p_motif_cur = p_motif_begin; p_motif_cur != NULL; p_motif_cur = p_motif_cur->next)
    {
        fprintf(fp_output, "%s\t%s\tTotal reads: %d\t Percentage of the total reads: %d%%\n\n", p_motif_cur->name, p_motif_cur->sequence, p_motif_cur->abundance, p_motif_cur->abundance*100/num_total_reads);
        fprintf(fp_output, "sequence\t\tlength\tabundance\tpercentage\n");
        for( p_read_cur = p_motif_cur->p_node_read; p_read_cur != NULL; p_read_cur = p_read_cur->next)
        {
            fprintf(fp_output, "%s\t\t%d\t%d\t%.2f%%\n", p_read_cur->sequence, p_read_cur->length, p_read_cur->abundance, (float)p_read_cur->abundance*100/p_motif_cur->abundance);
        }
        fprintf(fp_output, "\n===============================================================\n");
    }
    
    fclose(fp_output);
    
	return 0;	
}

/***********************************************************************************************/

/* GENERAL IO FUNCTIONS */

char*  read_seq(char* p_new, FILE* fp_input)
{
	char *pos;
    
    if (fgetc(fp_input) == EOF) return NULL;
    
    while (fgetc(fp_input) != '\n');     //skip the first line
    
    if (fgets(p_new, LEN_SEQ, fp_input) == NULL) return NULL;
    if ((pos=strchr(p_new, '\n')) != NULL)           *pos = '\0';
	
    while (fgetc(fp_input) != '\n');     //skip the third line
    while (fgetc(fp_input) != '\n');     //skip the fourth line
    
	return p_new;
}

struct node_motif*  read_motif(FILE* fp_motif)
{
    struct  node_motif  *p_new;
    char*   pos;
    char    ch;
    
	if( (p_new = malloc(sizeof(struct node_motif))) == NULL)
    {
        fprintf(stderr, "Not enough memory!\n");
        exit (1);
    }
    
    if ((ch = fgetc(fp_motif)) == EOF) return NULL;
    
    while ( ch != '>')
    {
        if ((ch = fgetc(fp_motif)) == EOF) return NULL;
    }
    
    fgets(p_new->name, LEN_MOTIF_NAME , fp_motif);
    if ((pos=strchr(p_new->name, '\n')) != NULL)           *pos = '\0';
	
    fgets(p_new->sequence, LEN_MOTIF_SEQ, fp_motif);
	if ((pos=strchr(p_new->sequence, '\n')) != NULL)       *pos = '\0';

    p_new->abundance = 0 ;
    p_new->p_node_read = NULL;
    p_new->next =   NULL;
    
    return p_new;
}

Bool        is_match(char* read, char* motif)
{
    if ( strstr_mis( read, motif, MIS_ALLOW) != NULL) return TRUE;
    else      return FALSE;
}


int strcmp_mis ( char* s1, char* s2, int n)     /* s1 is the pointer to the subject string, s2 substring, n is the numbr of mismatches allowed */
                                                /* for return, -1 - no match, mismatches are more than n,  m>=0 - remaining mismatches allowrance */
                                                /* s2 need a \0 at the the end, no such requirment on s1 */
{
    while ( *s2 )
    {
        if ((*s1++ - *s2++) && !(n--)) break;
    }
    return n;
}

char*   strstr_mis (char* s1, char* s2, int n )     /* Returns a pointer to the first occurrence of str2 in str1, or a null pointer if str2 is not part of str1. */
                                                    /* n is the mismatch allowrance */
                                                    /* strstr_mis search the whole s1, and report the FIRST matching of s2 with the LOWEST mismatch number */

{
    char * pos = NULL;
    int min_mis = -1, mismatch;
    
    while (*s1)
    {
        mismatch = strcmp_mis(s1, s2, n);
        if (mismatch > min_mis) min_mis = mismatch, pos = s1;
        s1++;
    }
    
    return pos;
}

struct node_read*        link_read(char* read, struct node_read* p_read_begin)
{
    struct node_read   *p_cur, *p_new, *p_temp;
    int length ;
    int compare;
    
    p_temp = NULL;
    length = strlen(read);
    
    
    for(p_cur = p_read_begin; p_cur != NULL;  p_cur = p_cur->next )
    {
        
        if (length < p_cur->length) break;
        if ( length == p_cur->length)
        {
            compare = strcmp(read, p_cur->sequence);
            if (compare < 0) break;
            else if (compare == 0) {p_cur->abundance++; return p_read_begin;}
        }
        p_temp = p_cur;
    }
    
    //need to creat a new nod, insert before p_cur
    
    if( (p_new = malloc(sizeof(struct node_read))) == NULL)
    {
        fprintf(stderr, "Not enough memory!\n");
        exit (1);
    }

    strcpy(p_new->sequence, read);
    p_new->length = length;
    p_new->abundance = 1;
    p_new->next = p_cur;
    
    if(p_temp == NULL) return p_new;
    else { p_temp->next = p_new; return p_read_begin;}
}







































