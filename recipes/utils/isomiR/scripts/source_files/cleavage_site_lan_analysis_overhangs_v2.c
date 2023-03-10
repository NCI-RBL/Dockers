
// This program read reference file in the format of fasta. Based on aligment info in SAM files, it convert the reads to sequence cover the cleavage sites.
// The length of the cleavage site will be two times of Windowns_radius. It start from the -R to R, with the site of mapping (first nuclotide of read) as 1.
// This program only measure the 5' end of aligned reads.

//July-24th-2016
//analyze the surrounding sequence of cleavage site!
//upstream 6 should be a C
//downsteam 4 should be a G



// July-26-2016

// This is a program try to proive a comprehensive analysis of the deep sequencing results of RNase III cleavage products

// **************************************************************************************************************************************************************************
//  cleavage_site_lan_analysis.exe reference input_SAM_file  length_option(all/number) reads_option_unique(yes/no) reads_option_threshold(number) output_file_option(on/off)
// **************************************************************************************************************************************************************************

// Total 6 parameters

// input_SAM file:              only the prefix! xxxx but not xxxx.sam

// length_option                all   -   reads of all length are considered in the analysis
//                              number (positive) : only reads of this particular length are considered in the analysis

// reads_option_unique          no  - all reads (including duplication) are counted. So reads abundance matters!
//                              yes   - only unique reads (no duplication) will be considered!

// reads_option_threshold       number (positive) :  only reads (peak) with abundance >= number will be considered in analysis.

// output_file_option:          on  -   output seuqnce around the cleavage site (xxxx.txt, xxxx_5P.txt, xxxx_3P.txt)
//                              off -   no output of sequence around the cleavage site


// July-26-2018 v2
// Try to modify this program. filter the output based on how many "ends" of reads supporting each cleavage site.
// Also try to provide a comprehansive way to analyze "hot spot"
// Output those "hot spot" sequences for sequence logo analysis.

// analysis_mode:               1: everything counts to support a cleavage site
//                              2: only 5' of read counts
//                              3: only sense-strand mapping ends (5' and 3') count
//                              4: only antisense-strand mapping ends (5' and 3') count
//                              3: only 5' end of sense-strand mapping ends count


//November 19th, 2018 cleavage_site_lan_analysis_overhangs.c
//This is basically chaning to a new program to analyze the overhangs
//the rough idea is this 1) only allow certain length (22mer e.g.) 2) fill in the table (array) as previously designed 3) checked the table and determine if any reads (+ or -) can find a parterner reads to form the duplex we expected.

//November 20th, 2018 v2
//Changing from only count the kind of reads (unique sequence) to count in their abundance.
//Also, now the threshold really mean something. It will be used to filter low abundant reads.

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define RADIUS              10                  // radius of motif around the cleavage
#define	MAX_REF_SIZE		2000                // max size of the refernce (both FF and MBP smaller than this value)
#define MAX_LINE_SIZE       256                 // max size of characters of one line of SAM file
#define	LEN_SEQ             150                 // max size of the reads

#define LENGTH_OF_INTEREST             22                 // the length of reads which will be analyzed. For the first version, 22mer was selected.


/*   DATA STRUCTRUE DEFINITITION HERE   */

char        reference[MAX_REF_SIZE];        // array holding the whole reference sequence
int         reference_sta[MAX_REF_SIZE][5] = {0};         // 0: conditional sum; 1: 5' +strand; 2: 3' +strand; 3: 5' -strand; 4: 3' -strand;

/*                                         FUNTION PROTOTYPE HERE                          */


/*   GLOABLE VIRABLE IS HERE   */


/************************************************************************************************************************************/


int main(int argc, char* argv[])		// Five arguments allowed
{
	
    FILE               *fp_input_sam;                               // only prefix of the file name
    FILE               *fp_reference;
    int                 threshold;                                  // any positive number
    int                 reference_length = 0;
    
    if (argc != 4)
	{
		fprintf(stderr, "\n cleavage_site_lan_analysis_overhangs_v2.exe reference input_SAM_file threshold(postive number)\n"); // the sam file name needs .sam at the end
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

    threshold = atoi(argv[3]);
    if ( threshold < 0 ) {fprintf(stderr, "threshold should be a positive number!\n"); exit (1);}        // error!

    
    //All the input arguments are correct!
    //*************************************************************************
    
    char    ch;
    int     index;
    
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
    }
    
    reference[index] = '\0';                                             // reference[index] is the reference sequence itself. It is 1_based! the first one reference[0] is N and the end is \0
    reference_length = index;
    
    //Reference sequence is loaded in reference[]
    //*************************************************************************
 
    
    char                    line[MAX_LINE_SIZE];           // one line of SAM file
    int                     flag, pos, length;
    char                    sequence[LEN_SEQ+1];           // hold the sequence of read
    int                     i;
    int                     num_cleavage_site = 0;

    
    while (fgets(line, MAX_LINE_SIZE, fp_input_sam) != NULL)
    {
        
        if (line[0] == '@') continue;
    
        sscanf(line, "%*s%d%*s%d%*s%*s%*s%*s%*s%s", &flag, &pos, sequence);
        
        length = strlen (sequence);
        
        if (length != LENGTH_OF_INTEREST) continue;         // only analyze 22mer
        
        if (flag == 0 )      // mapped to the +strand
        {
            reference_sta[pos][1]++;                //5' end
            reference_sta[pos+length][2]++;         //3' end
        }
        else if (flag == 16)    // mapped to the -strand
        {
            reference_sta[pos+length+2][3]++;       //5' end
            reference_sta[pos+2][4]++;              //3' end
        }
    }
    
   
    
    int     num_plus=0;
    int     num_minus=0;
    int     num_plus_yes=0;
    int     num_minus_yes=0;
    int     num_total=0;
    
    
    
    for(i=1; i<reference_length;i++)
    {
        reference_sta[i][0] = reference_sta[i][1] + reference_sta[i][2] + reference_sta[i][3] +reference_sta[i][4];
        
        if (reference_sta[i][1] != 0 && reference_sta[i][1] >= threshold)       // more than threshold..
        {
            num_plus+= reference_sta[i][1];
            if (reference_sta[i][4] != 0) num_plus_yes+= reference_sta[i][1];               // this position is suported by both plus(5) and minus strand(3')
        }
    }
    
    for(i=1; i<reference_length;i++)
    {
        if (reference_sta[i][3] != 0 && reference_sta[i][3] >= threshold)
        {
            num_minus+= reference_sta[i][3];
            if (reference_sta[i][2] != 0) num_minus_yes+= reference_sta[i][3];               // this position is suported by both plus(5) and minus strand(3')
        }
    }
    
    num_total = num_plus + num_minus;
    
    fprintf(stdout, "Pos\tCount\t5+\t3+\t5-\t3-\n");
    for(i=1; i<reference_length;i++)
    {
        if (reference_sta[i][0] >= threshold)
        {
            num_cleavage_site++;
            printf("%d\t%d\t%d\t%d\t%d\t%d\n",i, reference_sta[i][0], reference_sta[i][1], reference_sta[i][2], reference_sta[i][3], reference_sta[i][4] );
        }
    }
    
    fprintf(stdout, "Number of the reads mapped to plus strand is %d, among which %d are paired, that is %.3f%%\n", num_plus, num_plus_yes, (float)num_plus_yes*100/(float)num_plus);
    fprintf(stdout, "Number of the reads mapped to minus strand is %d, among which %d are paired, that is %.3f%%\n", num_minus, num_minus_yes, (float)num_minus_yes*100/(float)num_minus);
    fprintf(stdout, "Total Number of reads mapped is %d, among which %d are paired, that is %.3f%%\n", num_total, (num_plus_yes+num_minus_yes), (float)(num_plus_yes+num_minus_yes)*100/(float)num_total);
    
    fclose(fp_reference);
    fclose(fp_input_sam);
    
    return 0;
    
}


    


















