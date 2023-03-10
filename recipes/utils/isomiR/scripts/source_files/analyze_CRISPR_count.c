// 2016 Nov 17th

// ***************************************************************************************************************************************
// *    analyze_CRISPR_count.exe ~/Destop/NGS/index/sgRNA_list.txt input.sam output.txt                                            *
// ***************************************************************************************************************************************

// based on sam result to generate gene count!

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAX_LINE_SIZE   256
#define MAX_GENENAME_SIZE   100

struct gene                     // use to hold the list of reference genes.
{
	char            name[MAX_GENENAME_SIZE];
	int             count;
	struct gene *   next;
};


int main(int argc, char* argv[])		// three arguments allowed
{
	
    FILE                *fp_gene_list, *fp_sam, *fp_output;
    struct gene*        p_new = NULL;
    struct gene*        p_pos = NULL;
    struct gene*        p_cur = NULL;
    struct gene*        p_start = NULL;
    struct gene*        p_last = NULL;

    struct gene*        hash_table[100] = {NULL};        // cover range from A-Z, then h.
    
    char                gene_name[MAX_GENENAME_SIZE];
    char                line[MAX_LINE_SIZE];
    int                 num_gene=0, num_sgRNA=0, num_unmatched=0;
    int                 flag, unfound, i;
    char                first_letter;
    char                first_letter_cur;
    
    
    
    if (argc != 4)                              // three arguments!
	{
		fprintf(stderr, " analyze_CRISPR_count.exe ~/Destop/NGS/index/sgRNA_list.txt input.sam output.txt  \n");
		exit (1);
	}
    
    if ( (fp_gene_list = fopen(argv[1], "r")) == NULL )
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[1]);
        exit (1);
    }

    if ( (fp_sam = fopen(argv[2], "r")) == NULL )
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[2]);
        exit (1);
    }
    
    fp_output = fopen(argv[3], "w");
    
    
    //reading in data from sgRNA_list ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    if( (p_start = malloc(sizeof(struct gene))) == NULL)
    {
        fprintf(stderr, "Not enough memory!\n");
        exit (1);
    }
    p_start->next = NULL;
    p_last = p_start;
    
    while (fgets(gene_name, MAX_GENENAME_SIZE, fp_gene_list) != NULL )
    {
        if( (p_new = malloc(sizeof(struct gene))) == NULL)
        {
            fprintf(stderr, "Not enough memory!\n");
            exit (1);
        }
        
        i=0;
        while (gene_name[i] != '\n' && gene_name[i] != 13 && i < MAX_GENENAME_SIZE) {i++;}
        gene_name[i] = '\0';
        
        strcpy ( p_new->name, gene_name);
        p_new->count = 0;
        p_new->next = NULL;
        
        p_last->next = p_new;
        p_last = p_new;
        
        num_gene++;
    }
    
    // now the node chain of gene_list is ready. Note: the start node is empty! the real start is p_start->next!
    
    // generating hash_table!
    
    first_letter_cur = '\0';
    
    for (p_cur = p_start->next; p_cur != NULL; p_cur = p_cur->next)
    {
        first_letter = p_cur->name[0];
        if ( first_letter != first_letter_cur )
        {
            hash_table[first_letter - 'A'] = p_cur;
            first_letter_cur = first_letter;
        }
    }
    
    // start reading the sam file
    
    
    while (fgets(line, MAX_LINE_SIZE, fp_sam) != NULL)
    {
        
        if (line[0] == '@') continue;
        
        sscanf(line, "%*s%d%s%*d%*s%*s%*s%*s%*s%*s", &flag, gene_name);
        
        if (flag == 0 )      // mapped reads. updating the gene_name_list node chain...
        {
            num_sgRNA++;
            unfound = 1;         // set flag.
            
            first_letter = gene_name[0];
            p_pos = hash_table[first_letter-'A'];
            for (p_cur = p_pos; p_cur != NULL; p_cur = p_cur->next)
            {
                if(strcmp(gene_name, p_cur->name) == 0) { p_cur->count++;unfound = 0;break;} //find match
            }
            if (unfound) {num_unmatched++; printf("first letter is %c, and gene_name is %s", first_letter, gene_name); printf("hashtable point to %s\n", p_pos->name);}
        }
    }

    // print out the results
    
    fprintf(fp_output, "\t\tnumber of sgRNA\t%d\n", num_gene);
    fprintf(fp_output, "\t\tnumber of reads\t%d\n", num_sgRNA);
    fprintf(fp_output, "\t\tnumber of sgRNA reads not being analyzed\t%d\n", num_unmatched);
    
    for (p_cur = p_start->next; p_cur != NULL; p_cur = p_cur->next)
    {
        fprintf(fp_output, "%s\t%d\n", p_cur->name, p_cur->count);
    }

    fprintf(stderr, "\t\ttotal number of sgRNA is: %d\n", num_gene);
    fprintf(stderr, "\t\ttotal number of reads is: %d\n", num_sgRNA);
    fprintf(stderr, "\t\ttotal number of sgRNA reads not being analyzed is: %d\n", num_unmatched);
    
    fclose(fp_gene_list);
    fclose(fp_sam);
    fclose(fp_output);
    
    return 0;
    
}


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
