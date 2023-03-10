// 2015

// ***************************************************************************************************************************************
// *    cleavage_site_lan.exe Window_radius(R) reference input_SAM_file outupt_file                                                         *
// ***************************************************************************************************************************************

// This program read reference file in the format of fasta. Based on aligment info in SAM files, it convert the reads to sequence cover the cleavage sites.
// The length of the cleavage site will be two times of Windowns_radius. It start from the -R to R, with the site of mapping (first nuclotide of read) as 1.
// This program only measure the 5' end of aligned reads.

//July-24th-2016
//analyze the surrounding sequence of cleavage site!
//upstream 6 should be a C
//downsteam 4 should be a G
// ***************************************************************************************************************************************
// *    cleavage_site_lan_new.exe reference input_SAM_file                                                         *
// ***************************************************************************************************************************************


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
	
    FILE               *fp_reference, *fp_input_sam;
    char                ch;
    char                reference[MAX_REF_SIZE];        // array holding the whole reference sequence
    char*               ref_F_index[MAX_REF_SIZE];        // array holding pointers to strings//
    char*               ref_R_index[MAX_REF_SIZE];
    int                 index, size=0;
    
    
    if (argc != 3)
	{
		fprintf(stderr, "\n cleavage_site_lan_new.exe reference input_SAM_file  \n");
		exit (1);
	}
    
    if ( (fp_reference = fopen(argv[1], "r")) == NULL )
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[1]);
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
    
    char line[MAX_LINE_SIZE];
    int  flag, pos, length;
    char read[MAX_READ_SIZE];
    
    if ( (fp_input_sam = fopen(argv[2], "r")) == NULL )
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[2]);
        exit (1);
    }
    
    
    int         total=0;
    int         S5=0, W5=0, N5=0, S3=0, W3=0, N3=0;
    int         read_5S3S=0, read_5S3W=0, read_5S3N=0;
    int         read_5W3S=0, read_5W3W=0, read_5W3N=0;
    int         read_5N3S=0, read_5N3W=0, read_5N3N=0;
    int         flag_5, flag_3;
    
    while (fgets(line, MAX_LINE_SIZE, fp_input_sam) != NULL)
    {
        
        if (line[0] == '@') continue;
    
        sscanf(line, "%*s%d%*s%d%*s%*s%*s%*s%*s%s", &flag, &pos, read);
        
        length = strlen (read);
        
        if (flag == 0)      // mapped to the + strand, the pos the start of 5' end
        {
            total++;
            
            if (reference[pos-6] == 'C' && reference[pos+3] == 'G') {flag_5 = 2; S5++;}  // pos-6 is 5' upstream and pos+3 is 5' downstream
            else if (reference[pos-6] == 'C' || reference[pos+3] == 'G') {flag_5 = 1; W5++;}
            else {flag_5 = 0; N5++;}
            
            if (reference[pos-6+length] == 'C' && reference[pos+3+length] == 'G') {flag_3 = 2; S3++;}  // pos-6+length is 3' upstream and pos+3+length is 3' downstream
            else if (reference[pos-6+length] == 'C' || reference[pos+3+length] == 'G') {flag_3 = 1; W3++;}
            else {flag_3 = 0; N3++;}
            
            if (flag_5 == 2)
            {
                if(flag_3 == 2) read_5S3S++;
                else if (flag_3 == 1) read_5S3W++;
                else read_5S3N++;
            }
            else if(flag_5 ==1)
            {
                if(flag_3 == 2) read_5W3S++;
                else if (flag_3 == 1) read_5W3W++;
                else read_5W3N++;
            }
            else
            {
                if(flag_3 == 2) read_5N3S++;
                else if (flag_3 == 1) read_5N3W++;
                else read_5N3N++;
            }
        }
        
        if (flag == 16)   //mapped to the - strand, the pos the start of 3'end
        {
            total++;
            
            if (reference[pos+5] == 'G' && reference[pos-4] == 'C') {flag_3 = 2; S3++;}  // pos+5 is 3' upstream and pos-4 is 3' downstream
            else if (reference[pos+5] == 'G' || reference[pos-4] == 'C') {flag_3 = 1; W3++;}
            else {flag_3 = 0; N3++;}
            
            if (reference[pos+5+length] == 'G' && reference[pos-4+length] == 'C') {flag_5 = 2; S5++;}  // pos+5+length is 5' upstream and pos-4+length is 5' downstream
            else if (reference[pos+5+length] == 'G' || reference[pos-4+length] == 'C') {flag_5 = 1; W5++;}
            else {flag_5 = 0; N5++;}
            
            if (flag_5 == 2)
            {
                if(flag_3 == 2) read_5S3S++;
                else if (flag_3 == 1) read_5S3W++;
                else read_5S3N++;
            }
            else if(flag_5 ==1)
            {
                if(flag_3 == 2) read_5W3S++;
                else if (flag_3 == 1) read_5W3W++;
                else read_5W3N++;
            }
            else
            {
                if(flag_3 == 2) read_5N3S++;
                else if (flag_3 == 1) read_5N3W++;
                else read_5N3N++;
            }

            
        }
    }
    
    //printf("total number of reads is: %d\n", total);
    
    //printf("5S\t5W\t5N\n");
    //printf("%4d%%\t%4d%%\t%4d%%\n", S5*100/total, W5*100/total, N5*100/total);

    //printf("3S\t3W\t3N\n");
    //printf("%4d%%\t%4d%%\t%4d%%\n", S3*100/total, W3*100/total, N3*100/total);
    
    //printf("5S3S\t5S3W\t5S3N\n");
    //printf("%4d%%\t%4d%%\t%4d%%\n", read_5S3S*100/total, read_5S3W*100/total, read_5S3N*100/total);
    //printf("5W3S\t5W3W\t5W3N\n");
    //printf("%4d%%\t%4d%%\t%4d%%\n", read_5W3S*100/total, read_5W3W*100/total, read_5W3N*100/total);
    //printf("5N3S\t5N3W\t5N3N\n");
    //printf("%4d%%\t%4d%%\t%4d%%\n", read_5N3S*100/total, read_5N3W*100/total, read_5N3N*100/total);
    
    printf("total, S5, W5, N5, S3, W3, N3, read_5S3S, read_5S3W, read_5S3N, read_5W3S, read_5W3W, read_5W3N, read_5N3S, read_5N3W, read_5N3N\n");
    printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", total, S5, W5, N5, S3, W3, N3, read_5S3S, read_5S3W, read_5S3N, read_5W3S, read_5W3W, read_5W3N, read_5N3S, read_5N3W, read_5N3N );
    printf("\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\n",  S5*100/total, W5*100/total, N5*100/total, S3*100/total, W3*100/total, N3*100/total, read_5S3S*100/total, read_5S3W*100/total, read_5S3N*100/total, read_5W3S*100/total, read_5W3W*100/total, read_5W3N*100/total, read_5N3S*100/total, read_5N3W*100/total, read_5N3N*100/total );
    
    fclose(fp_reference);
    fclose(fp_input_sam);
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    