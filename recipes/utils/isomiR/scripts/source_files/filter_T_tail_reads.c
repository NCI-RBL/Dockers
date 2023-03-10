// 2017 September 19

//  This program take ready reads (xxx_ready.fastq) as input. Sorting out reads contain T tails and some statistics..

// ***************************************************************************************************************************************
// *    filter_T_tail_reads.exe Min_number_T pass_rate input.fastq output_contain_T.fastq output_no_T_tail.fastq                         *
// ***************************************************************************************************************************************


// allow input pass_rate as the an argument -- type float (between 70~100)


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define Bool		int
#define	TRUE		1
#define	FALSE       0
#define SUCCESS     1
#define FAIL        0


#define	LEN_SEQ		150


#define PASS_RATE       100           // used in searching motif in reads, no mismatch allowed. this is the default, will be overwritten by command line argument


/*   DATA STRUCTRUE DEFINITITION HERE   */

struct fastq
{
    char        name[LEN_SEQ];
    char        sequence[LEN_SEQ];
    char		tag[LEN_SEQ];
    char		score[LEN_SEQ];
};

/****************************************/


/*                                         FUNTION PROTOTYPE HERE                          */

Bool    read_fastq(struct fastq *, FILE*);
void    print_fastq(struct fastq *, FILE*);

Bool    Has_T_tail(struct fastq *, int min_number_T, int pass_rate);


/************************************************************************************************************************************/


int main(int argc, char* argv[])		// Five arguments allowed
{
	
    FILE                *fp_input, *fp_output_yes, *fp_output_no;
    long    int         num_total_reads = 0;
    long    int         num_Ttail_yes = 0;
    long    int         num_Ttail_no = 0;
    int                 min_number_T = 0;
    int                 pass_rate = PASS_RATE;
    
    struct fastq        *p_read;
    struct fastq        read;
    
    if (argc != 6)
	{
		fprintf(stderr, "\n filter_T_tail_reads.exe Min_number_T pass_rate input.fastq output_contain_T.fastq output_no_T_tail.fastq \n");
		exit (1);
	}
    
    min_number_T = atoi(argv[1]);
    if ( min_number_T < 0 || min_number_T >= 50)
    {
        fprintf(stderr, " Invalid minimal number of T in tail, should be a number betweeen 0 and 50\n");
        exit (1);
    }
    
    pass_rate = atoi(argv[2]);
    if ( pass_rate < 50 || pass_rate > 100)
    {
        fprintf(stderr, " Invalid pass rate, should be a number betweeen 50 and 100\n");
        exit (1);
    }
    
    if ( (fp_input = fopen(argv[3], "r")) == NULL )                     // open motif file, making sure it is sucessful..
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[3]);
        exit (1);
    }

	fp_output_yes = fopen(argv[4], "w");
    fp_output_no  = fopen(argv[5], "w");
    
    p_read = & read;
    
    while ( read_fastq(p_read, fp_input) != FAIL)
    {
        num_total_reads++;
        
        if (Has_T_tail(p_read, min_number_T, pass_rate))
        {
            print_fastq(p_read, fp_output_yes);
            num_Ttail_yes++;
        }
        else
        {
            print_fastq(p_read, fp_output_no);
            num_Ttail_no++;
        }
    }
    
    fclose(fp_input);
    fclose(fp_output_yes);
    fclose(fp_output_no);
    
    printf("Number of total reads is %ld\n",num_total_reads);
    printf("Number of reads having a Tail with minimal %d T and a passing rate of %d%% is %ld, that is %.2f%% of total reads\n", min_number_T, pass_rate, num_Ttail_yes, (float)num_Ttail_yes*100/(float)num_total_reads);
	printf("Number of reads without such tail is %ld, that is %ld%% of total reads\n", num_Ttail_no, num_Ttail_no*100/num_total_reads);
    
    
	return 0;	
}

/***********************************************************************************************/

/* GENERAL IO FUNCTIONS */



Bool    Has_T_tail(struct fastq * p_read, int min_number_T, int pass_rate)
{
    int length = strlen(p_read->sequence);
    int number_T=0;
    int number = 0;
    float correct_rate;
    char    *p;
    
    p = p_read->sequence + length;      // p point to the last '\0'
    p = p -2;                           // p point to the last nucleotide
    
    //printf ("the last one is %c\n", *p);
    
    while ( p != p_read->sequence - 1)
    {
        number ++;
        if (*p == 'T') number_T++;
        
        correct_rate = ((float) number_T)/((float) number) * 100;        // it is OK for read only containing part of the motif, as long as it over the pass_rate
    
        if (correct_rate >= pass_rate && number_T >= min_number_T) return TRUE;
        
        p--;
    }

    return FALSE;
}


Bool    read_fastq(struct fastq* p_read, FILE* fp_input)
{
    if (fgets(p_read->name, LEN_SEQ, fp_input) == NULL) return FAIL;
    if (fgets(p_read->sequence, LEN_SEQ, fp_input) == NULL) return FAIL;
    if (fgets(p_read->tag, LEN_SEQ, fp_input) == NULL) return FAIL;
    if (fgets(p_read->score, LEN_SEQ, fp_input) == NULL) return FAIL;
    
    return SUCCESS;
}

void print_fastq(struct fastq *p_read, FILE* fp_output)
{
    fprintf(fp_output, "%s", p_read->name);
    fprintf(fp_output, "%s", p_read->sequence);
    fprintf(fp_output, "%s", p_read->tag);
    fprintf(fp_output, "%s", p_read->score);
}







































