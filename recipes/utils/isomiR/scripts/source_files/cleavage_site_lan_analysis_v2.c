
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






#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define RADIUS              10                  // radius of motif around the cleavage
#define	MAX_REF_SIZE		2000                // max size of the refernce (both FF and MBP smaller than this value)
#define MAX_LINE_SIZE       256                 // max size of characters of one line of SAM file
#define	LEN_SEQ             150                 // max size of the reads




/*   DATA STRUCTRUE DEFINITITION HERE   */

char        reference[MAX_REF_SIZE];        // array holding the whole reference sequence
int         reference_sta[MAX_REF_SIZE][5] = {0};         // 0: conditional sum; 1: 5' +strand; 2: 3' +strand; 3: 5' -strand; 4: 3' -strand;

/*                                         FUNTION PROTOTYPE HERE                          */

int              print_cleavage_motif(int pos);

/*   GLOABLE VIRABLE IS HERE   */

FILE        *fp_output      = NULL;

/************************************************************************************************************************************/


int main(int argc, char* argv[])		// Five arguments allowed
{
	
    FILE               *fp_input_sam;                               // only prefix of the file name
    FILE               *fp_reference;
    char                file_name[256];
    char                sam_file_name[256];
    char                output_file_name[256];
    int                 output_file_option;                         // on: 1 off:0
    int                 mode;                                       // 1, 2, 3, 4 or 5.
    int                 threshold;                                  // any positive number
    int                 reference_length = 0;
    
    if (argc != 6)
	{
		fprintf(stderr, "\n cleavage_site_lan_analysis_overhangs.exe reference input_SAM_file threshold(postive number)\n");
		exit (1);
	}
    
    if ( (fp_reference = fopen(argv[1], "r")) == NULL )
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[1]);
        exit (1);
    }

    
    strcpy (file_name, argv[2]);
    strcpy (sam_file_name, file_name);
    strcat (sam_file_name, ".sam");
    if ( (fp_input_sam = fopen(sam_file_name, "r")) == NULL )
    {
        fprintf(stderr, "Cannot find file %s!\n", sam_file_name);
        exit (1);
    }

    mode = atoi(argv[3]);
    if (mode < 1 || mode > 5)
    {fprintf(stderr, "analysis_mode should be between 1 and 5!\n"); exit (1);}        // error!
    
    
    threshold = atoi(argv[4]);
    if ( threshold < 0 ) {fprintf(stderr, "threshold should be a positive number!\n"); exit (1);}        // error!
    
   
    if ( strcmp (argv[5], "off") == 0 ) output_file_option = 0;                                  // off is 0
    else if ( strcmp (argv[5], "on") == 0 )
    {
        output_file_option = 1;                                                                 // on is 1
        
        strcpy(output_file_name, file_name);
        strcat(output_file_name, ".txt");
        
        fp_output = fopen(output_file_name, "w");

    }
    else {fprintf(stderr, "output_file_option should be on or off!\n"); exit (1);}        // error!
    
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
    int                     skipped = 0;
    
    while (fgets(line, MAX_LINE_SIZE, fp_input_sam) != NULL)
    {
        
        if (line[0] == '@') continue;
    
        sscanf(line, "%*s%d%*s%d%*s%*s%*s%*s%*s%s", &flag, &pos, sequence);
        
        length = strlen (sequence);
        
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
    
    for(i=1; i<reference_length;i++)
    {
        if (mode == 1) reference_sta[i][0] = reference_sta[i][1] + reference_sta[i][2] + reference_sta[i][3] +reference_sta[i][4];
        else if (mode == 2) reference_sta[i][0] = reference_sta[i][1] + reference_sta[i][3];
        else if (mode == 3) reference_sta[i][0] = reference_sta[i][1] + reference_sta[i][2];
        else if (mode == 4) reference_sta[i][0] = reference_sta[i][3] + reference_sta[i][4];
        else if (mode == 5) reference_sta[i][0] = reference_sta[i][1];
        else {fprintf(stderr, "Illigal mode number!\n"); exit (1);}
    }
    
    fprintf(stdout, "Pos\tCount\t5+\t3+\t5-\t3-\n");
    for(i=1; i<reference_length;i++)
    {
        if (reference_sta[i][0] >= threshold)
        {
            num_cleavage_site++;
            printf("%d\t%d\t%d\t%d\t%d\t%d\n",i, reference_sta[i][0], reference_sta[i][1], reference_sta[i][2], reference_sta[i][3], reference_sta[i][4] );
            if (output_file_option == 1) {if (print_cleavage_motif(i) == 1) skipped++;};
        }
    }
    
    fprintf(stdout, "Number of the cleavage sites is %d\n", num_cleavage_site);
    fprintf(stdout, "Number of the cleavage sites too close to the ends is %d\n", skipped);
    
    if (output_file_option == 1)
    {
        fclose(fp_output);
    }
    fclose(fp_reference);
    fclose(fp_input_sam);
    
    return 0;
    
}

int    print_cleavage_motif(int pos)                               // return 0 if sucessful, return 1 if too close to the ends!
{
    char motif[2*RADIUS+1];
    
    int ref_size;
    ref_size = strlen (reference);
    
    if ( pos <= RADIUS || pos + RADIUS >= ref_size) return 1;    // the pos is too close to the end of reference!
        
    strncpy(motif, reference+pos-RADIUS, 2*RADIUS );           // motif from -R to R (no 0), 2R long. and pos point to the +1 position.
    motif[2*RADIUS] = '\0';
    
    fprintf(fp_output, "%s\n",motif);

    return 0;
}

    


















