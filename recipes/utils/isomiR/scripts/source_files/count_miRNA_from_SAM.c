// Aug-29-2018

// This is a program try to count the number of miRNA reads based on SAM mapping to the "hairpin"!

// **************************************************************************************************************************************************************************
//  count_miRNA_from_SAM.exe reference input_SAM_file
// **************************************************************************************************************************************************************************

//      1. read in refernce file, build a node tree to hold it.
//      2. calculate the middel point position of hairpin, which will be used later to infer a mature miRNA sequence is from 5p or from 3p.
//      3. read in SAM file, chekcing the flag of mapping
//      4. if "0" -> mapped, determine if 5p or 3p, then revise the reference table
//      5. other value -> skip
//      6. report the note tree.


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define LEN_HAIRPIN         250                 // max size of the miRNA hairpin
#define LEN_MIRNA_NAME      60                 // max size of the miRNA name
#define MAX_LINE_SIZE       256                 // max size of characters of one line of SAM file
#define LEN_SEQ             150                 // max length of the miRNA reads

#define RADIUS              10                  // radius of motif around the cleavage
#define	MAX_REF_SIZE		2000                // max size of the refernce (both FF and MBP smaller than this value)




/*   DATA STRUCTRUE DEFINITITION HERE   */

struct miRNA
{
    char        miRNA_name[LEN_MIRNA_NAME+1];
    char        miRNA_sequence[LEN_HAIRPIN +1];
    int         middle;
    int         num_5p;
    int         num_3p;
    struct miRNA *     next;
};

/*  FUNTION PROTOTYPE HERE  */

/*   GLOABLE VIRABLE IS HERE   */


/************************************************************************************************************************************/


int main(int argc, char* argv[])		// Five arguments allowed
{
	
    FILE               *fp_reference;
    FILE               *fp_input_sam;
    
    if (argc != 3)
	{
        fprintf(stderr, "\n count_miRNA_from_SAM.exe reference input_SAM_file\n");
		exit (1);
	}
    
    if ( (fp_reference = fopen(argv[1], "r")) == NULL )
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[1]);
        exit (1);
    }
    
    if ( (fp_input_sam = fopen(argv[2], "r")) == NULL )
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[2]);
        exit (1);
    }
    
    //All the input arguments are correct!
    //*************************************************************************
    
    char                    line[MAX_LINE_SIZE];           // one line of the fasta file
    struct miRNA            *p_miRNA_start = NULL;
    struct miRNA            *p_miRNA_new;
    struct miRNA            *p_miRNA_cur;
    
    int                     i;
    int                     num_miRNA;
    
    while (fgets(line, MAX_LINE_SIZE, fp_reference) != NULL)
    {
        if( (p_miRNA_new = malloc(sizeof(struct miRNA))) == NULL)
        {
            fprintf(stderr, "Not enough memory!\n");
            exit (1);
        }
        
    
        if (line[0] != '>') {fprintf(stderr, "reading fa file missing the '>' !\n"); exit (1);}
            
        for (i=1; i < LEN_MIRNA_NAME; i++)
            {
                if ( line[i] == ' ') break;
                else p_miRNA_new->miRNA_name[i-1] = line[i];
            }
        p_miRNA_new->miRNA_name[i-1] = '\0';                                // Note: No '\n' at the end!
        
        if (fgets(line, MAX_LINE_SIZE, fp_reference) == NULL)
        {
            fprintf(stderr, "Missing the second line in the fasta file!\n");
            exit (1);
        }
       
        strcpy(p_miRNA_new->miRNA_sequence, line);
        p_miRNA_new->middle = strlen(p_miRNA_new->miRNA_sequence)/2;
        
        p_miRNA_new->num_5p = 0;
        p_miRNA_new->num_3p = 0;
        p_miRNA_new->next = NULL;
        
        if (p_miRNA_start != NULL)
        {
            p_miRNA_cur->next = p_miRNA_new;
        }
        else
        {
            p_miRNA_start = p_miRNA_new;
        }
        
        p_miRNA_cur = p_miRNA_new;
        num_miRNA++;
    }
    
    //Reference sequence is loaded in reference[]
    //*************************************************************************
    
    int                     flag, pos, length;
    char                    sequence[LEN_SEQ+1];           // hold the sequence of read
    char                    miRNA_name_ref[LEN_MIRNA_NAME+1];       // hold the name of miRNA the reads mapped to
    
 
    
    while (fgets(line, MAX_LINE_SIZE, fp_input_sam) != NULL)
    {
        
        if (line[0] == '@') continue;
    
        sscanf(line, "%*s%d%s%d%*s%*s%*s%*s%*s%s", &flag,  miRNA_name_ref, &pos, sequence);
        
        length = strlen (sequence);
        
        if (flag != 0 ) continue;       // unmapped reads, skip to read next line in SAM file
        
        for (p_miRNA_cur = p_miRNA_start; p_miRNA_cur != NULL; p_miRNA_cur = p_miRNA_cur->next)
        {
            if (strcmp(p_miRNA_cur->miRNA_name,miRNA_name_ref) == 0)
            {
                if ( pos < p_miRNA_cur->middle) {p_miRNA_cur->num_5p++;}
                else {p_miRNA_cur->num_3p++;}
                break;
            }
        }
    }
    
    for (p_miRNA_cur = p_miRNA_start; p_miRNA_cur != NULL; p_miRNA_cur = p_miRNA_cur->next)
    {
        fprintf(stdout, "%s-5p", p_miRNA_cur->miRNA_name);
        fprintf(stdout, "\t%d\n", p_miRNA_cur->num_5p);
        fprintf(stdout, "%s-3p", p_miRNA_cur->miRNA_name);
        fprintf(stdout, "\t%d\n", p_miRNA_cur->num_3p);
    }
    
    fprintf(stdout, "Number of the miRNA in reference file is %d\n", num_miRNA);
    
 
    fclose(fp_reference);
    fclose(fp_input_sam);
    
    return 0;
    
}













