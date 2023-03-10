// July-25-2018

// The purpose of this program is to help analyze RNase III cleavage reads for Lan's project.
// The program take fastq format file (matched reads) as input and generate the following outputs:
//      1.  Reads aligned at the 5' end (txt format), the length is defined as OUTPUT_LENGTH (15 in this version). reads shorter than this will have N filled at the 3'end
//      2.  Reeds aligned at the 3' end, reads shorter than the OUTPUT_LENGTH will have N filled at the 5' end
//      3.  Statistics of these reads, including total number of reads, number of unique reads of both input and output.
//      4.  Statistics of the aligments. How many A, T, G, C at each position.
//      5.  Filters including:  a) read abundance (taken as input argument) and b) output collapsed reads (unique reads) or all reads (read abundance matters!)
// Of note, this program could be modified in the future to collapse reads.

// **************************************************************************************************************************************************************************
//  cleavage_reads_analyze_lan.exe  input_FASTQ_file  output_reads_option_threshold(number) output_reads_option_unique(yes/no)  output_file_option(on/off)
// **************************************************************************************************************************************************************************

// Total 4 parameters

// input_FASTQ file:                    only the prefix! xxxx but not xxxx.fastq

// output_reads_option_threshold        number (positive) :  only reads (peak) with abundance >= number will be considered in output and analysis.

// reads_option_unique                  no  - all reads (including duplication) are counted. So reads abundance matters!
//                                      yes   - only unique reads (no duplication) will be included in the outputs!

// output_file_option:                  on  -   output sequence alignments ( xxxx_5P.txt, xxxx_3P.txt)
//                                      off -   no output of sequence alignments


// July-14th-2016
// This is for Lan to find motif in all cleave sites
// Purpose of this program is to collapse reads and then output reads higher than a given number. In addition, it will some statistics regarding the distribution of read abundance.
// the input is in the format of txt -> one read per line. need to change code to read fastq or fasta format in the future. the output is also the txt format.
// " collapse_reads.exe min_number_reads_to_be_reported (INT) input output

 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define 	Bool		int
#define     TRUE		1
#define     FALSE       0

#define	LEN_NAME	150
#define	LEN_SEQ		150
#define	LEN_TAG		25

#define OUTPUT_LENGTH              10                  // length of sequence in output alignments



/*   DATA STRUCTRUE DEFINITITION HERE   */

struct node_read
{
	char                sequence[LEN_SEQ+1];
	int					abundance; 
	struct node_read * 	next;
};

struct fastq
{
    char     name[LEN_NAME+1];
    char     sequence[LEN_SEQ+1];
    char    tag[LEN_TAG+1];
    char    score[LEN_SEQ+1];
};

/*  FUNTION PROTOTYPE HERE   */

struct node_read * add_node ( char * );  /* return pointer to the newly added node */

char*    read_fastq_sequence(char* read);

void analyze_reads (struct node_read* p_node_cur);


/*   GLOABLE VIRABLE IS HERE   */

FILE        *fp_input_fastq    = NULL;
FILE        *fp_output_5p      = NULL;
FILE        *fp_output_3p      = NULL;

int                     output_file_option;                                // on: 1 off:0
int                     output_reads_option_unique;                        // yes:1 no:0
int                     output_reads_option_threshold;                     // yes:positive number no:0

int         num_input_total=0;                  // total number of input reads
int         num_input_unique=0;                 // total number of unique input reads
int         num_output_total=0;                  // total number of output reads
int         num_output_unique=0;                 // total number of unique output reads

int         alignment_5p[OUTPUT_LENGTH][4];             // hold 5p alignment statistic information, 0-A, 1-T, 2-G, 3-C
int         alignment_3p[OUTPUT_LENGTH][4];             // hold 3p alignment statistic information, 0-A, 1-T, 2-G, 3-C


int main(int argc, char* argv[])		/* Four arguments */
{
    char                    file_name[256];
    char                    fastq_file_name[256];
    char                    output_file_name[256];
    char                    output_file_name_5p[256];
    char                    output_file_name_3p[256];
    
    
    if (argc != 5)
	{
		fprintf(stderr, "\n cleavage_reads_analyze_lan.exe  input_FASTQ_file  output_reads_option_threshold(number) output_reads_option_unique(yes/no)  output_file_option(on/off) \n");
		exit (1);
	}
    
    strcpy (file_name, argv[1]);
    strcpy (fastq_file_name, file_name);
    strcat (fastq_file_name, ".fastq");
    if ( (fp_input_fastq = fopen(fastq_file_name, "r")) == NULL )
    {
        fprintf(stderr, "Cannot find file %s!\n", fastq_file_name);
        exit (1);
    }
    
    output_reads_option_threshold = atoi(argv[2]);
    if ( output_reads_option_threshold < 0 ) {fprintf(stderr, "output_reads_option_threshold should be a positive number!\n"); exit (1);}        // error!
    
   
    if (strcmp(argv[3], "yes") == 0) output_reads_option_unique = 1;
    else if (strcmp(argv[3], "no") == 0) output_reads_option_unique = 0;
    else {fprintf(stderr, "output_reads_option_unique should be yes or no!\n"); exit (1);}
    
    if ( strcmp (argv[4], "off") == 0 ) output_file_option = 0;                                  // off is 0
    else if ( strcmp (argv[4], "on") == 0 )
    {
        output_file_option = 1;                                                                 // on is 1
        
        strcpy(output_file_name, file_name);
        if (output_reads_option_unique) strcat(output_file_name, "_unique");
        strcat(output_file_name, "_equal_or_more_than_");
        strcat(output_file_name, argv[2]);
        
        strcpy(output_file_name_5p, output_file_name);
        strcpy(output_file_name_3p, output_file_name);
        strcat(output_file_name_5p, "_5P_alignment.txt");
        strcat(output_file_name_3p, "_3P_alignment.txt");
        
        fp_output_5p = fopen(output_file_name_5p, "w");
        fp_output_3p = fopen(output_file_name_3p, "w");
    }
    else {fprintf(stderr, "output_file_option should be on or off!\n"); exit (1);}        // error!
    
    //All the input arguments are correct!
    //*************************************************************************
    
    char                    read[LEN_SEQ+1];             // hold each read sequence
    
    struct node_read        *p_node_start = NULL;        // p_node_start point to the first node
    struct node_read        *p_node_cur = NULL;
    

    
    while (read_fastq_sequence(read) != NULL)
    {
        num_input_total++;
        
      //  fprintf(stdout, "sequence of the current input reads: %s", read);
        
        if ( p_node_start == NULL) p_node_start = add_node (read);
        else
        {
            for (p_node_cur = p_node_start; p_node_cur != NULL; p_node_cur = p_node_cur->next)
            {
                if(strcmp(read, p_node_cur->sequence) == 0) { p_node_cur->abundance++;break;} //find match
                
                if(p_node_cur->next == NULL) { p_node_cur->next = add_node(read);break;}      // reach the end of node tree
            }
        }
        
    }
    
    // at this point, all reads have been stored in the node tree
    
    for ( p_node_cur = p_node_start; p_node_cur != NULL; p_node_cur = p_node_cur->next)
    {
        num_input_unique++;
        if (p_node_cur->abundance >= output_reads_option_threshold) analyze_reads (p_node_cur);
    }

    fprintf(stdout, "\nTotal number of input reads:\t%d\n", num_input_total);
    fprintf(stdout, "Total number of unique input reads:\t%d\n", num_input_unique);
    fprintf(stdout, "\n");
    fprintf(stdout, "reads of a minimal abundance of %d are analyzed\n", output_reads_option_threshold);
    if (output_reads_option_unique == 1) fprintf(stdout, "outputs are unique reads only\n");
    fprintf(stdout, "Total number of output reads:\t%d\n", num_output_total);
    fprintf(stdout, "Total number of unique output reads:\t%d\n", num_output_unique);
    fprintf(stdout, "\n");
    
	fprintf(stdout, "Reads are aligned at the 5' end\n");
    fprintf(stdout, "Pos\tA\tT\tG\tC\n");
    
    int i;
    for (   i = 0; i<OUTPUT_LENGTH; i++)
    {
        fprintf(stdout, "%d\t%d\t%d\t%d\t%d\n", i, alignment_5p[i][0], alignment_5p[i][1],alignment_5p[i][2],alignment_5p[i][3]);
    }
    
    fprintf(stdout, "\n");
    fprintf(stdout, "Reads are aligned at the 3' end\n");
    fprintf(stdout, "Pos\tA\tT\tG\tC\n");
    
    for (   i = 0; i<OUTPUT_LENGTH; i++)
    {
        fprintf(stdout, "%d\t%d\t%d\t%d\t%d\n", i, alignment_3p[i][0], alignment_3p[i][1],alignment_3p[i][2],alignment_3p[i][3]);
    }
    
    fclose (fp_input_fastq);
    fclose (fp_output_5p);
    fclose (fp_output_3p);
	return 0;	
}


void analyze_reads (struct node_read* p_node_cur)
{
    int     count = 0;
    if (output_reads_option_unique == 1) count = 1;
    else count = p_node_cur->abundance;
    
    num_output_unique++;
    num_output_total += count;
    
    char    motif_5p[OUTPUT_LENGTH+2];        //string holding the reads from 5p end, the last two chars will be \n then \0
    char    motif_3p[OUTPUT_LENGTH+2];
    
    int i;
    int flag = 0 ;
    
    for (i=0; i < OUTPUT_LENGTH; i++)
    {
        if (flag) {motif_5p[i] = 'N'; continue;}                                                        // read is shorter than the OUTPUT_LENGTH
        if      (p_node_cur->sequence[i] == 'A') { alignment_5p[i][0] += count; motif_5p[i] = 'A'; }
        else if (p_node_cur->sequence[i] == 'T') { alignment_5p[i][1] += count; motif_5p[i] = 'T'; }
        else if (p_node_cur->sequence[i] == 'G') { alignment_5p[i][2] += count; motif_5p[i] = 'G'; }
        else if (p_node_cur->sequence[i] == 'C') { alignment_5p[i][3] += count; motif_5p[i] = 'C'; }
        else if (p_node_cur->sequence[i] == '\n' || p_node_cur->sequence[i] == '\0')                    // read is shorter than the OUTPUT_LENGTH, reach the end, set flag
        {motif_5p[i] = 'N'; flag = 1;}
        else {fprintf(stderr, "Illigal character %c! Its code is %d\n", p_node_cur->sequence[i], p_node_cur->sequence[i]); exit (1);}        // error!
    }
    
    motif_5p[i] = '\n';
    motif_5p[i+1] = '\0';
    
    if (output_file_option)
    {
        for(i=0; i<count; i++)
        {fprintf(fp_output_5p, "%s", motif_5p);}
    }
    
    // now deal with 3p analysis and alignment
    
    int p ;
    int length;
    length = strlen(p_node_cur->sequence);
    p = length -2;                  // sequence[p] will be the last Char not including \n
    
    for (i=0; i < OUTPUT_LENGTH; i++)
    {
        p = (length - 2) - i;                               // read is shorter than the OUTPUT_LENGTH
        if (p < 0) {motif_3p[i] = 'N';}                     // reach the begining of the read!
        else if (p_node_cur->sequence[p] == 'A') { alignment_3p[i][0] += count; motif_3p[i] = 'A'; }
        else if (p_node_cur->sequence[p] == 'T') { alignment_3p[i][1] += count; motif_3p[i] = 'T'; }
        else if (p_node_cur->sequence[p] == 'G') { alignment_3p[i][2] += count; motif_3p[i] = 'G'; }
        else if (p_node_cur->sequence[p] == 'C') { alignment_3p[i][3] += count; motif_3p[i] = 'C'; }
        else {fprintf(stderr, "Illigal character %c! Its code is %d\n", p_node_cur->sequence[i], p_node_cur->sequence[i]); exit (1);}        // error!
    }
    
    motif_3p[i] = '\n';
    motif_3p[i+1] = '\0';
    
    if (output_file_option)
    {
        for(i=0; i<count; i++)
        {fprintf(fp_output_3p, "%s", motif_3p);}
    }
    
}

char*    read_fastq_sequence(char* read)
{
    struct fastq temp;
    
    if (fgets(temp.name, LEN_NAME, fp_input_fastq ) == NULL) return NULL;
    if (fgets(temp.sequence, LEN_SEQ, fp_input_fastq) == NULL) return NULL;
    if (fgets(temp.tag, LEN_TAG, fp_input_fastq) == NULL) return NULL;
    if (fgets(temp.score, LEN_SEQ, fp_input_fastq) == NULL) return NULL;
    
    strcpy (read, temp.sequence);
/*
    fprintf(stdout, "\nname:  %s", temp.name);
    fprintf(stdout, "sequence:  %s", temp.sequence);
    fprintf(stdout, "tag:  %s", temp.tag);
    fprintf(stdout, "score:  %s", temp.score);
*/
    return read;
}

struct node_read * add_node ( char * read)
{
	struct node_read *p_new;
	
	if( (p_new = malloc(sizeof(struct node_read))) == NULL)
		{
			fprintf(stderr, "Not enough memory!\n");
			exit (1);
		}
	
	strcpy (p_new->sequence, read);
	p_new->abundance = 1;
	p_new->next = NULL;
    
	return p_new;	
}




