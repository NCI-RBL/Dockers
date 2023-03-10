// 2017 September 19

//  This program take ready reads (xxx_ready.fastq) as input. Sorting out reads contain T tails and some statistics..

// ***************************************************************************************************************************************
// *    filter_T_tail_reads_v2.exe Min_number_T pass_rate input.fastq output_contain_T.fastq output_no_T_tail.fastq                         *
// ***************************************************************************************************************************************


// allow input pass_rate as the an argument -- type float (between 70~100)

// 2017 September 20  v2
// compare two fastq files (WT sv DIS3KO)
// make min_number_T and pass_rate variable...


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
float   T_tail_rate( FILE* , int , int );

/************************************************************************************************************************************/


int main(int argc, char* argv[])		// Five arguments allowed
{
	
    FILE                *fp_input_WT, *fp_input_KO;
    int                 min_number_T = 0;
    int                 pass_rate = PASS_RATE;
    
    struct fastq        *p_read;
    struct fastq        read;
    
    if (argc != 3)
	{
		fprintf(stderr, "\n filter_T_tail_reads.exe input_WT.fastq input_KO.fastq \n");
		exit (1);
	}
    
    
    if ( (fp_input_WT = fopen(argv[1], "r")) == NULL )                     // open motif file, making sure it is sucessful..
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[1]);
        exit (1);
    }
    
    if ( (fp_input_KO = fopen(argv[2], "r")) == NULL )                     // open motif file, making sure it is sucessful..
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[2]);
        exit (1);
    }
    
    fclose(fp_input_WT);
    fclose(fp_input_KO);

    
    printf ("%s\n", argv[1]);
    printf ("Pass_rate\t");
    for (pass_rate = 75; pass_rate <=100; pass_rate++)
    {
        printf ("%d\t", pass_rate);
    }
    printf ("\n");
    
    for (min_number_T = 3; min_number_T < 20; min_number_T++)
    {
        printf ("%d\t", min_number_T);
        
        for (pass_rate = 75; pass_rate <=100; pass_rate++)
        {
            fp_input_WT = fopen(argv[1], "r");
            printf ("%.6f\t", T_tail_rate(fp_input_WT, min_number_T, pass_rate));
            fclose(fp_input_WT);
        }
        printf ("\n");
    }
    printf ("\n");
    printf ("\n");
    
    printf ("%s\n", argv[2]);
    printf ("Pass_rate\t");
    for (pass_rate = 75; pass_rate <=100; pass_rate++)
    {
        printf ("%d\t", pass_rate);
    }
    printf ("\n");
    
    for (min_number_T = 3; min_number_T < 20; min_number_T++)
    {
        printf ("%d\t", min_number_T);
        
        for (pass_rate = 75; pass_rate <=100; pass_rate++)
        {
            fp_input_KO = fopen(argv[2], "r");
            printf ("%.6f\t", T_tail_rate(fp_input_KO, min_number_T, pass_rate));
            fclose(fp_input_KO);
        }
        printf ("\n");
    }
    
    
    
	return 0;	
}

/***********************************************************************************************/

/* GENERAL IO FUNCTIONS */

float   T_tail_rate( FILE* fp_input, int min_number_T, int pass_rate)
{
    long    int         num_total_reads = 0;
    long    int         num_Ttail_yes = 0;
    
    struct fastq        *p_read;
    struct fastq        read;

    p_read = & read;
    
    while ( read_fastq(p_read, fp_input) != FAIL)
    {
        num_total_reads++;
        if (Has_T_tail(p_read, min_number_T, pass_rate)) num_Ttail_yes++;
    }
    
    return  (float)num_Ttail_yes*100/(float)num_total_reads;
}


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







































