// 2015-Aug-31st :

// ***************************************************************************************************************************************
// *    cleavage_site_lan.exe Window_radius(R) reference input_SAM_file outupt_file                                                         *
// ***************************************************************************************************************************************

// This program read reference file in the format of fasta. Based on aligment info in SAM files, it convert the reads to sequence cover the cleavage sites.
// The length of the cleavage site will be two times of Windowns_radius. It start from the -R to R, with the site of mapping (first nuclotide of read) as 1.
// This program only measure the 5' end of aligned reads.


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define Bool		int
#define	TRUE		1
#define	FALSE       0
#define SUCCESS     1
#define FAIL        0


#define	MAX_REF_SIZE		2000
#define MAX_LINE_SIZE       255
#define MAX_READ_SIZE       100




/*   DATA STRUCTRUE DEFINITITION HERE   */


/****************************************/


/*                                         FUNTION PROTOTYPE HERE                          */

void rev_com(char* d, char* s, int length);

/************************************************************************************************************************************/


int main(int argc, char* argv[])		// Four arguments allowed
{
	
    FILE               *fp_reference, *fp_input_sam, *fp_output;
    int                 r, length;
    char                ch;
    char                reference[MAX_REF_SIZE];        // array holding the whole reference sequence
    char*               ref_F_index[MAX_REF_SIZE];        // array holding pointers to strings//
    char*               ref_R_index[MAX_REF_SIZE];
    int                 index, size=0;
    
    
    if (argc != 5)
	{
		fprintf(stderr, "\n cleavage_site_lan.exe Window_radius(R) reference input_SAM_file outupt_file \n");
		exit (1);
	}
    
    r = atoi(argv[1]);
    if ( r <= 0 || r > 20)
    {
        fprintf(stderr, " Invalid radius of windows, should be a number betweeen 0 and 20\n");
        exit (1);
    }
    
    
    if ( (fp_reference = fopen(argv[2], "r")) == NULL )                     // open motif file, making sure it is sucessful..
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[2]);
        exit (1);
    }

    
    while  ((ch = fgetc(fp_reference)) != '\n')                             // skp the first line >....
    {
        ;
    }
    
    reference[0] = 'N';
    
    for (index = 1; (ch = fgetc(fp_reference)) != EOF; index++)
    {
        if (index == MAX_REF_SIZE)
        {
            
            fprintf(stderr, "\n refernce is too large! \n");
            exit (1);
        }
        reference[index] = ch;
        size++;
    }
    
    reference[index] = '\0';                                             // reference[index] is the reference sequence itself. It is 1_based! the first one reference[0] is N and the end is \0
    
    length = 2*r + 2;                                                   // length-1 is the size of output motif
    
    for (index = r + 1; index <= size - r; index++)
    {
        ref_F_index[index] = malloc(length);
        ref_R_index[index] = malloc(length);
        
        strncpy(ref_F_index[index], reference+index-r, length -1);
        ref_F_index[index][length-1] = '\0';
        
        rev_com(ref_R_index[index], ref_F_index[index], length);
    }
    
//    for (index = r + 1; index <=size -r; index++ )
//    {
//        printf("%s\n", ref_F_index[index]);
//        printf("%s\n", ref_R_index[index]);
//    }
    
    char line[MAX_LINE_SIZE];
    int  flag, pos, read_length;
    char read[MAX_READ_SIZE];
    
    if ( (fp_input_sam = fopen(argv[3], "r")) == NULL )                     // open motif file, making sure it is sucessful..
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[3]);
        exit (1);
    }

    fp_output = fopen(argv[4], "w");
  
  //  fprintf(fp_output, "%s\n", reference);
  //  fprintf(fp_output, "size is %d\n", size);
  //  fprintf(fp_output, "r is %d\n", r);
  //  fprintf(fp_output, "length is %d\n", length);
    
    while (fgets(line, MAX_LINE_SIZE, fp_input_sam) != NULL)
    {
        
        if (line[0] == '@') continue;
    
        sscanf(line, "%*s%d%*s%d%*s%*s%*s%*s%*s%s", &flag, &pos, read);
        
        if (flag == 0)
        {
            if (pos >= (r + 1) && pos <= (size - r)) fprintf(fp_output, "%s\n",ref_F_index[pos]);
        }
        
        if (flag == 16)
        {
            read_length = strlen(read);
            pos = pos + read_length - 1;
            if (pos >= (r + 1) && pos <= (size - r)) fprintf(fp_output, "%s\n",ref_R_index[pos]);
        }
    }
    
    fclose(fp_reference);
    fclose(fp_input_sam);
    fclose(fp_output);
    
    return 0;
    
}


void rev_com(char* d, char* s, int length)
{
    int i;
    
    for (i = 0; i < length - 1; i++)
    {
        if (s[length - 2 - i ] == 'A' ) { d[i] = 'T'; continue;}
        if (s[length - 2 - i ] == 'T' ) { d[i] = 'A'; continue;}
        if (s[length - 2 - i ] == 'G' ) { d[i] = 'C'; continue;}
        if (s[length - 2 - i ] == 'C' ) { d[i] = 'G'; continue;}
    }
}
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    