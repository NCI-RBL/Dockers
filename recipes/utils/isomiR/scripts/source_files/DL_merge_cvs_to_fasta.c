// 2016 Nov 15th

// ***************************************************************************************************************************************
// *    DL_merge_cvs_to_fasta.exe input_1.cvs input_2.cvs output.fa                                                                      *
// ***************************************************************************************************************************************

// used to merge two sgRNA libaray to fasta format, which can then be used as refrence file for bowtie..

// 2016 Nov 18th

// add two more outputs. 1. list of all sgRNA names 2. combined gene_name_list
// ***************************************************************************************************************************************
// *    DL_merge_cvs_to_fasta.exe input_1.cvs input_2.cvs output.fa output_sgRNA_name.txt output_gene_name.txt                           *
// ***************************************************************************************************************************************


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAX_LINE_SIZE   100

struct gene
{
	char            ID[50];
	char            seq[25];
	int             index;
	struct gene *   next;
};

struct gene *       read_gene(FILE*);

int main(int argc, char* argv[])		// three arguments allowed
{
	
    FILE                *fp_input_A, *fp_input_B, *fp_output_fa, *fp_output_sgRNA, *fp_output_gene;
    struct gene*        p_cur = NULL;
    struct gene*        p_start = NULL;
    struct gene*        p_last = NULL;
    struct gene*        p = NULL;
    struct gene*        p_position = NULL;
    struct gene*        p_temp = NULL;
    
    char                ch;
    int                 num_gene=0, num_sgRNA=0, num_gene_printed=0;
    char                gene_name[50] ={'\0'};
    
    if (argc != 6)
	{
		fprintf(stderr, "\n DL_merge_cvs_to_fasta.exe input_1.cvs input_2.cvs output.fa  output_sgRNA_name.txt output_gene_name.txt\n");
		exit (1);
	}
    
    if ( (fp_input_A = fopen(argv[1], "r")) == NULL )
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[1]);
        exit (1);
    }

    if ( (fp_input_B = fopen(argv[2], "r")) == NULL )
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[2]);
        exit (1);
    }
    
    fp_output_fa = fopen(argv[3], "w");
    fp_output_sgRNA = fopen(argv[4], "w");
    fp_output_gene = fopen(argv[5], "w");
    
    while  ((ch = fgetc(fp_input_A)) != '\n')                             // skp the first line for file1
    {
        ;
    }
    
    while  ((ch = fgetc(fp_input_B)) != '\n')                             // skp the first line for file2
    {
        ;
    }
    
    
    //reading in data from library_A ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    if ((p_start = read_gene(fp_input_A)) == NULL) exit (1);
    p_start->index = 1;
    num_gene = 1;
    p_last=p_start;
    
    while ((p_cur = read_gene(fp_input_A)) != NULL )
    {
        p_last->next = p_cur;                                                           // added to the node chain
        
        if ( strcmp (p_cur->ID, p_last->ID) == 0)   p_cur->index = p_last->index+1;          // same gene..
        else                                        {p_cur->index = 1; num_gene++;}                      // different gene..
        p_last = p_cur;
    }
    
    printf("done reading A\n");
    //reading in data from library_B ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    p_position = p_start;
    
    while ((p_cur = read_gene(fp_input_B)) != NULL )
    {
        
        for(p = p_position; p != NULL; p = p->next)
        {
            if ( strcmp ( p_cur->ID, p->ID ) == 0 )
            {
                if (p->next == NULL) break;
                if ( strcmp (p_cur->ID, (p->next)->ID) != 0)
                {
                    p_temp = p->next;
                    p->next = p_cur;
                    p_cur->next = p_temp;
                    
                    p_cur->index = p->index +1;
                    p_position = p_cur;
                    break;
                }
            }
        }
    }

    printf("done reading B\n");
    
    for (p = p_start; p != NULL; p = p->next)
    {
        
        if (p->ID[1] != 'o')
        {
            fprintf(fp_output_fa, ">%s_%d\n%s\n", p->ID, p->index, p->seq);
            fprintf(fp_output_sgRNA,"%s_%d\n", p->ID, p->index);
        }
        else                            // Non-target controls...
        {
            fprintf(fp_output_fa, ">%s\n%s\n", p->ID, p->seq);
            fprintf(fp_output_sgRNA,"%s\n", p->ID);
        }
        if( strcmp (gene_name, p->ID) != 0)
        {
            fprintf(fp_output_gene,"%s\n", p->ID);
            strcpy ( gene_name, p->ID);
            num_gene_printed++;
        }
        
        
        num_sgRNA++;
    }
    
        

    
    fprintf(stderr, "total number of sgRNA is: %d\n", num_sgRNA);
    fprintf(stderr, "total number of gene is: %d\n", num_gene);
    fprintf(stderr, "total number of gene printed is: %d\n", num_gene_printed);
    
    fclose(fp_input_A);
    fclose(fp_input_B);
    fclose(fp_output_fa);
    fclose(fp_output_sgRNA);
    fclose(fp_output_gene);
    
    return 0;
    
}

struct gene *       read_gene(FILE*fp)
{
	struct              gene  *p_new;
    char                line[MAX_LINE_SIZE];                 // array holding one line in input files
    
    if( (p_new = malloc(sizeof(struct gene))) == NULL)
    {
        fprintf(stderr, "Not enough memory!\n");
        exit (1);
    }
    
    if (fgets(line, MAX_LINE_SIZE, fp) == NULL) return NULL;
    
    sscanf(line, "%s%*s%s", p_new->ID, p_new->seq);
    
    p_new->next = NULL;
    
    return p_new;
}

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    