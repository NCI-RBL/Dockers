// 2015-Mar-30th : v2
// file name changed to analyzed_reads_by_motif_v2
// add in two more argument in command line: The first is minimal percentage for a read to be reported; the second will be the total percentage to be reported.
// only when BOTH conditions fail, the report will not be reported!
// both arguments are float numbers
// Changed the way we measure mismatch -- using pass_rate
// allow input pass_rate as the third arguments -- type float (between 70~100)

// ***************************************************************************************************************************************
// *    analyze_read_by_motif_v2.exe Min_percentage_to_be_reported Max_percentage_in_total pass_rate motif_list input.fastq output.txt   *
// ***************************************************************************************************************************************

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

#define PASS_RATE       100           // used in searching motif in reads, no mismatch allowed. this is the default, will be overwritten by command line argument


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


struct node_motif*  read_motif(FILE* fp_motif);             //read one motif (two lines) from the file fp_motif pointing to. return a pointer if sucessful, otherwise return a NULL pointer.

char*               read_seq(char* read_cur, FILE* fp_input);               //read one sequence from fastq file.. return NULL when reaching end...

Bool        is_match(char* read, char* motif, float pass_rate);             // using pass_rate to determine how many mismatches are allowed..

struct node_read*        link_read(char* read, struct node_read* p_read_begin);  //link the read to the node_read tree

void   sort_and_print( FILE* fp_output, struct node_read* p_read_begin, float min_percentage, float max_percentage, float total); // sort by abundance, print to output file when condistions met

/************************************************************************************************************************************/


int main(int argc, char* argv[])		// Six arguments allowed
{
	
    FILE               *fp_motif, *fp_input, *fp_output;
    int                 num_total_reads = 0;
    struct node_motif*  p_motif_begin = NULL, *p_motif_cur;
    struct node_read*   p_read_cur;
    float               min_percentage, max_percentage, pass_rate = PASS_RATE;
    
    if (argc != 7)
	{
		fprintf(stderr, "\n analyze_read_by_motif_v2.exe Min_percentage_to_be_reported Max_percentage_in_total pass_rate motif_list input.fastq output.txt \n");
		exit (1);
	}
    
    min_percentage = (float) (atof(argv[1]));
    if ( min_percentage < 0 || min_percentage >= 50)
    {
        fprintf(stderr, " Invalid minimal percentage, should be a number betweeen 0 and 50\n");
        exit (1);
    }
    
    max_percentage = (float) (atof(argv[2]));
    if ( max_percentage < 50 || max_percentage > 100)
    {
        fprintf(stderr, " Invalid maximal total percentage, should be a number betweeen 50 and 100\n");
        exit (1);
    }
    
    pass_rate = (float) (atof(argv[3]));
    if ( pass_rate < 70 || pass_rate > 100)
    {
        fprintf(stderr, " Invalid passing rate, should be a number betweeen 70 and 100\n");
        exit (1);
    }
    
    
    if ( (fp_motif = fopen(argv[4], "r")) == NULL )                     // open motif file, making sure it is sucessful..
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
        fprintf(stderr, "\n %s file is empty!\n", argv[4]);
        exit (1);
    }

    if ( (fp_input = fopen(argv[5], "r")) == NULL )                     // open input file, making sure it is sucessful..
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[5]);
        exit (1);
    }

    char read_cur[LEN_SEQ + 1] = {'\0'};                        // maximal number of char read in is LEN_SEQ, So the last '\0' is always there..
    
    while ( read_seq(read_cur, fp_input) != NULL)
    {
        num_total_reads++;
        
        for( p_motif_cur = p_motif_begin; p_motif_cur != NULL; p_motif_cur = p_motif_cur->next)
        {
          
            if ( is_match( read_cur, p_motif_cur->sequence, pass_rate))
            {
                p_motif_cur->abundance++;
                p_motif_cur->p_node_read = link_read(read_cur, p_motif_cur->p_node_read);
                break;
            }
        }
    }
    
    fclose(fp_input);
    
    fp_output = fopen(argv[6], "w");
    
    fprintf(fp_output, "Motif list file is %s; Input file is %s; Number of total reads is %d\n", argv[4], argv[5], num_total_reads);
    fprintf(fp_output, "Minimal percentage to be reported is %s%%; Minimal total percentage is %s%%; the passing rate is %s%%!\n\n", argv[1], argv[2],argv[3]);
	
    for( p_motif_cur = p_motif_begin; p_motif_cur != NULL; p_motif_cur = p_motif_cur->next)
    {
        fprintf(fp_output, "%s\t%s\tTotal reads: %d\t Percentage of the total reads: %d%%\n\n", p_motif_cur->name, p_motif_cur->sequence, p_motif_cur->abundance, p_motif_cur->abundance*100/num_total_reads);
        fprintf(fp_output, "sequence\t\tlength\tabundance\tpercentage\n");
     
        sort_and_print(fp_output, p_motif_cur->p_node_read, min_percentage, max_percentage, (float) p_motif_cur->abundance);
        
        fprintf(fp_output, "\n===============================================================\n");
    }
    
    fclose(fp_output);
    
	return 0;	
}

/***********************************************************************************************/

/* GENERAL IO FUNCTIONS */

void   sort_and_print( FILE* fp_output, struct node_read* p_read_begin, float min_percentage, float max_percentage, float total)
{
    
    struct node_read *p_read_cur, *p_read_prev, *p_max, *p_max_prev;
    int max;
    float percentage, total_percentage = 0;
    
    
    while (p_read_begin != NULL)
    {
        for( p_read_cur = p_read_begin, p_read_prev = NULL, max = 0; p_read_cur != NULL; p_read_prev = p_read_cur, p_read_cur = p_read_cur->next)          // try to find the read with max abundance
        {
            if (p_read_cur->abundance > max)
            {
                max = p_read_cur->abundance;
                p_max = p_read_cur;
                p_max_prev = p_read_prev;
            }
        }
    
        percentage = (float)p_max->abundance*100/total;
        if (total_percentage >= max_percentage && percentage < min_percentage)
            {
                fprintf(fp_output, "Total percentage\t\t\t\t\t%.2f%%\n", total_percentage);
                return;
            }
        fprintf(fp_output, "%s\t\t%d\t%d\t%.2f%%\n", p_max->sequence, p_max->length, p_max->abundance, percentage);
        total_percentage += percentage;
    
        if (p_max_prev == NULL) p_read_begin = p_max->next;
        else    p_max_prev->next = p_max->next;
    }
    
    fprintf(fp_output, "Total percentage\t\t\t\t\t%.2f%%\n", total_percentage);
    
    return;
}


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

Bool        is_match(char* read, char* motif, float pass_rate)
{
    int length = strlen(motif);
    int match;
    float correct_rate;
    char    *p_read, *p_motif;
    
    while ( *read )
    {
        p_read = read;
        p_motif = motif;
        match = 0;
        while (*p_read && *p_motif)
        {
            if (*p_read == *p_motif) match++;
            p_read++;
            p_motif++;
        }
        
        correct_rate = ((float) match)/((float) length) * 100;        // it is OK for read only containing part of the motif, as long as it over the pass_rate
    
        if (correct_rate >= pass_rate ) return TRUE;
        
        read++;
    }
    
    return FALSE;
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







































