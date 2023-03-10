// 2016 Nov 20th

// ***************************************************************************************************************************************
// *    analyze_by_gene.exe input_counts_by_sgRNA.txt output_counts_by_gene.txt                                          *
// ***************************************************************************************************************************************

// convert counts by sgRNA to counts by gene.

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAX_LINE_SIZE   256
#define MAX_GENENAME_SIZE   50



int main(int argc, char* argv[])		// three arguments allowed
{
	
    FILE                *fp_sgRNA, *fp_output;
    
    
    char                line[MAX_LINE_SIZE];
    char                gene_name[MAX_GENENAME_SIZE];
    char                gene_name_cur[MAX_GENENAME_SIZE] = {'\0'};
    
    int                 count, total= -1;
    int                 num_gene=0, num_read=0;
    
    char*               pos;
    
    
    
    
    if (argc != 3)                              // three arguments!
	{
		fprintf(stderr, " analyze_by_gene.exe input_counts_by_sgRNA.txt output_counts_by_gene.txt\n");
		exit (1);
	}
    
    if ( (fp_sgRNA= fopen(argv[1], "r")) == NULL )
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[1]);
        exit (1);
    }

    fp_output = fopen(argv[2], "w");
    
    
    //reading in data from counts by sgRNA ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    while (fgets(line, MAX_LINE_SIZE, fp_sgRNA) != NULL)
    {
        
        if (line[0] == '\t') continue;
        
        sscanf(line, "%s%d", gene_name, &count);
        
        pos = gene_name;
        
        while ( *pos != '_' && *pos != '\0') {pos++;}
        *pos = '\0';
        
        if ( strcmp (gene_name, gene_name_cur) != 0)
        {
            if (total >= 0) {fprintf(fp_output, "\t%d\n", total); num_read += total; }       // if not the first time, then print the total count for the last gene.
            fprintf(fp_output, "%s", gene_name);
            strcpy(gene_name_cur, gene_name);
            total = count;
            
            num_gene++;
        }
        else total += count;
    }
    
    fprintf(fp_output, "\t%d\n", total);  // print out the last gene's count.
    num_read += total;
        
        
    
    fprintf(fp_output, "\t\tnumber of gene\t%d\n", num_gene);
    fprintf(fp_output, "\t\tnumber of reads\t%d\n", num_read);
    
    fprintf(stderr, "\t\tnumber of gene\t%d\n", num_gene);
    fprintf(stderr, "\t\tnumber of reads\t%d\n", num_read);
    
    fclose(fp_sgRNA);
    fclose(fp_output);
    
    return 0;
    
}


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    