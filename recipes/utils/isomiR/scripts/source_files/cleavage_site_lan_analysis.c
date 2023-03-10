
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


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define RADIUS              10                  // radius of motif around the cleavage
#define	MAX_REF_SIZE		2000                // max size of the refernce (both FF and MBP smaller than this value)
#define MAX_LINE_SIZE       256                 // max size of characters of one line of SAM file
#define	LEN_SEQ             75                  // max size of the reads


/*   DATA STRUCTRUE DEFINITITION HERE   */

struct read_sam                                 // contain info of reads in SAM file
{
	int                 Flag;
    int                 Pos;
    char                Sequence[LEN_SEQ+1];
    int                 Sequence_length;
    int                 Abundance;
    struct read_sam     *next;
};


/*                                         FUNTION PROTOTYPE HERE                          */

struct  read_sam *  add_node (int   flag, int pos, char* sequence, int read_length);
void                print_cleavage_motif(struct read_sam * p, int output_file_option);
void                analyze_read(struct read_sam * p);
void                rev_com(char* s);            // reverse compliment the string!

/*   GLOABLE VIRABLE IS HERE   */

FILE        *fp_output      = NULL;
FILE        *fp_output_5p      = NULL;
FILE        *fp_output_3p      = NULL;

char        reference[MAX_REF_SIZE];        // array holding the whole reference sequence
int         ref_size = 0;

int         total=0;                // total number of reads mapped to the reference
int         num_unique_peak =0;         // total number of unique reads
int         num_unique_peak_after_length = 0;
int         num_unique_peak_analyzed =0;
int         reads_analyzed =0;         // total number of reads being analyzed


int         S5=0, W5=0, N5=0, S3=0, W3=0, N3=0;
int         read_5S3S=0, read_5S3W=0, read_5S3N=0;
int         read_5W3S=0, read_5W3W=0, read_5W3N=0;
int         read_5N3S=0, read_5N3W=0, read_5N3N=0;
int         flag_5, flag_3;
int         G5=0, G3=0, C5=0, C3=0;

/************************************************************************************************************************************/


int main(int argc, char* argv[])		// Six arguments allowed
{
	
    FILE               *fp_input_sam;                               // only prefix of the file name
    FILE               *fp_reference;
    char                file_name[256];
    char                sam_file_name[256];
    char                output_file_name[256];
    char                output_file_name_5p[256];
    char                output_file_name_3p[256];
    int                 output_file_option;                         // on: 1 off:0
    int                 length_option;                              // 0 - all, otherwise, reads with this specific length will be analyzed
    int                 reads_option_unique;                        // yes:1 no:0
    int                 reads_option_threshold;                     // yes:1 no:0
    
    if (argc != 7)
	{
		fprintf(stderr, "\n cleavage_site_lan_analysis.exe reference input_SAM_file  length_option(all/number) reads_option_unique(yes/no) reads_option_threshold(number) output_file_option(on/off)\n");
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

    if (strcmp(argv[3], "all") == 0) length_option = 0;             // all length -> length_option = 0
    else
    {
            
        length_option = atoi(argv[3]);
        if ( length_option <= 0 || length_option >= 50) {fprintf(stderr, "length_option should be between 0 and 50!\n"); exit (1);}        // error!
    }
    
    if (strcmp(argv[4], "yes") == 0) reads_option_unique = 1;
    else if (strcmp(argv[4], "no") == 0) reads_option_unique = 0;
    else {fprintf(stderr, "reads_option_unique should be yes or no!\n"); exit (1);}
    
    
    reads_option_threshold = atoi(argv[5]);
    if ( reads_option_threshold < 0 ) {fprintf(stderr, "reads_option_threshold should be a positive number!\n"); exit (1);}        // error!
    
    
    
    if ( strcmp (argv[6], "off") == 0 ) output_file_option = 0;                                  // off is 0
    else if ( strcmp (argv[6], "on") == 0 )
    {
        output_file_option = 1;                                                                 // on is 1
        
        strcpy(output_file_name, file_name);
        strcat(output_file_name, "_length_");
        strcat(output_file_name, argv[3]);
        if (reads_option_unique) strcat(output_file_name, "_unique");
        strcat(output_file_name, "_bigger_than_");
        strcat(output_file_name, argv[5]);
        
        strcpy(output_file_name_5p, output_file_name);
        strcpy(output_file_name_3p, output_file_name);
        strcat(output_file_name, ".txt");
        strcat(output_file_name_5p, "_cleavage_5P.txt");
        strcat(output_file_name_3p, "_cleavage_3P.txt");
        
        fp_output = fopen(output_file_name, "w");
        fp_output_5p = fopen(output_file_name_5p, "w");
        fp_output_3p = fopen(output_file_name_3p, "w");
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
        ref_size++;
    }
    
    reference[index] = '\0';                                             // reference[index] is the reference sequence itself. It is 1_based! the first one reference[0] is N and the end is \0
    
    
    //Reference sequence is loaded in reference[]
    //*************************************************************************
 
    
    char                    line[MAX_LINE_SIZE];           // one line of SAM file
    int                     flag, pos, read_length;
    char                    sequence[LEN_SEQ+1];           // hold the sequence of read
    struct read_sam         *p_node_start = NULL;		/* p_node_start point to the first node */
	struct read_sam         *p_node_cur = NULL;
    int                     i;
    
    while (fgets(line, MAX_LINE_SIZE, fp_input_sam) != NULL)
    {
        
        if (line[0] == '@') continue;
    
        sscanf(line, "%*s%d%*s%d%*s%*s%*s%*s%*s%s", &flag, &pos, sequence);
        
        read_length = strlen (sequence);
        
        if (flag == 0 || flag == 16)      // mapped reads. input into the node chain which will be analyized..
        {
            total++;
            if ( p_node_start == NULL) p_node_start = add_node (flag, pos, sequence, read_length);
            else
            {
                for (p_node_cur = p_node_start; p_node_cur != NULL; p_node_cur = p_node_cur->next)
                {
                    if(strcmp(sequence, p_node_cur->Sequence) == 0) { p_node_cur->Abundance++;break;} //find match
                    
                    if(p_node_cur->next == NULL) { p_node_cur->next = add_node (flag, pos, sequence, read_length);break;}      // reach the end of node tree
                }
            }
        }
    }
    
// at this point, all reads have been stored in the node tree. The tree is pointed by p_nod_start
//********************************************************************
    
    
    for (p_node_cur = p_node_start; p_node_cur != NULL; p_node_cur = p_node_cur->next)
    {
        num_unique_peak++;
        
        if(length_option != 0 && length_option != p_node_cur->Sequence_length) continue;          // reads cannot pass length_option requirment
        
        num_unique_peak_after_length++;
        
        if(p_node_cur->Abundance < reads_option_threshold ) continue;                    // reads without enough abundance (>= threshold) be to considered
        
        num_unique_peak_analyzed++;
        
        if (!reads_option_unique)
        {
            
            for (i = 0; i < p_node_cur->Abundance; i++)
            {
                print_cleavage_motif(p_node_cur, output_file_option);
                analyze_read(p_node_cur);
            }
        }
    
        if ( reads_option_unique )
        {
            print_cleavage_motif(p_node_cur, output_file_option);
            analyze_read(p_node_cur);
        }
    
    }
    
    if (length_option == 0)
    {
        printf("Unique\tThreshold\tLength\ttotal\tnum_unique_peak\treads_analyzed\tnum_unique_peak_after_length\tnum_unique_peak_analyzed\t\tS\tW\tG\tC\tN\t\tS5\tW5\tG5\tC5\tN5\t\tS3\tW3\tG3\tC3\tN3\t\t\tread_5S3S\tread_5S3W\tread_5S3N\tread_5W3S\tread_5W3W\tread_5W3N\tread_5N3S\tread_5N3W\tread_5N3N\n");
    
    
        printf("%s\t%s\t%s\t", argv[4], argv[5],argv[3]);
    
        printf("%d\t%d\t%d\t%d\t%d", total,num_unique_peak,  reads_analyzed, num_unique_peak_after_length, num_unique_peak_analyzed);
        printf("\t\t%d\t%d\t%d\t%d\t%d\t\t%d\t%d\t%d\t%d\t%d\t\t%d\t%d\t%d\t%d\t%d\t\t\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", (S5+S3)/2, (W5+W3)/2, (G5+G3)/2, (C5+C3)/2, (N5+N3)/2, S5, W5, G5, C5, N5, S3, W3, G3, C3, N3, read_5S3S, read_5S3W, read_5S3N, read_5W3S, read_5W3W, read_5W3N, read_5N3S, read_5N3W, read_5N3N );
    }
    
    printf("\t\t%s\t\t\t", argv[3]);
    
    if (reads_analyzed != 0)
    {
        
        printf("%d%%\t%d%%\t%d%%\t", reads_analyzed*100/total, num_unique_peak_after_length*100/num_unique_peak, num_unique_peak_analyzed*100/num_unique_peak );
        printf("\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t\t\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\n",  (S5+S3)*50/reads_analyzed, (W5+W3)*50/reads_analyzed, (G5+G3)*50/reads_analyzed, (C5+C3)*50/reads_analyzed,(N5+N3)*50/reads_analyzed, S5*100/reads_analyzed, W5*100/reads_analyzed, G5*100/reads_analyzed,C5*100/reads_analyzed, N5*100/reads_analyzed, S3*100/reads_analyzed, W3*100/reads_analyzed, G3*100/reads_analyzed,C3*100/reads_analyzed, N3*100/reads_analyzed, read_5S3S*100/reads_analyzed, read_5S3W*100/reads_analyzed, read_5S3N*100/reads_analyzed, read_5W3S*100/reads_analyzed, read_5W3W*100/reads_analyzed, read_5W3N*100/reads_analyzed, read_5N3S*100/reads_analyzed, read_5N3W*100/reads_analyzed, read_5N3N*100/reads_analyzed );
    }
    else printf("\n");
    
    fclose(fp_reference);
    fclose(fp_input_sam);
    
    if (output_file_option == 1)
    {
        fclose(fp_output);
        fclose(fp_output_5p);
        fclose(fp_output_3p);
    }
        
    return 0;
    
}


struct read_sam * add_node(int flag, int  pos,char* sequence, int read_length)
{
	struct read_sam *p_new;
	
	if( (p_new = malloc(sizeof(struct read_sam))) == NULL)
    {
        fprintf(stderr, "Not enough memory!\n");
        exit (1);
    }
	
	p_new->Flag = flag;
    p_new->Pos = pos;
    strcpy (p_new->Sequence, sequence);
	p_new->Sequence_length = read_length;
    p_new->Abundance = 1;
	p_new->next = NULL;
    
	return p_new;	
}

void    print_cleavage_motif(struct read_sam * p, int option)
{
    int pos, length, flag;
    char motif_5p[2*RADIUS+1];
    char motif_3p[2*RADIUS+1];
    
    pos = p->Pos;
    length = p->Sequence_length;
    flag = p->Flag;
    
    if (option == 0) return;
   
    if (flag == 16) pos = pos + length ;                 // mapped to the - strand, now pos point to the cleavage site of 5'end of this read (on positive strand)
        
    if ( pos <= RADIUS || pos + RADIUS >= ref_size) return;    // the pos is too close to the end of reference!
        
    strncpy(motif_5p, reference+pos-RADIUS, 2*RADIUS );           // motif from -R to R (no 0), 2R long. and pos point to the +1 position.
    motif_5p[2*RADIUS] = '\0';
        
    if (flag == 16) rev_com ( motif_5p);


    pos = p->Pos;                                        // reset pos to the original
    if (flag == 0) pos = pos + length ;                // mapped to the - strand, now pos point to the cleavage site of 5'end of this read (on positive strand)
            
    if ( pos <= RADIUS || pos + RADIUS >= ref_size) return;    // the pos is too close to the end of reference!
        
    strncpy(motif_3p, reference+pos-RADIUS, 2*RADIUS );           // motif from -R to R (no 0), 2R long. and pos point to the +1 position.
    motif_3p[2*RADIUS] = '\0';
        
    if (flag == 16) rev_com ( motif_3p);
    
    
 
    fprintf(fp_output, "%s\n",motif_5p);
    fprintf(fp_output, "%s\n",motif_3p);
    fprintf(fp_output_5p, "%s\n",motif_5p);
    fprintf(fp_output_3p, "%s\n",motif_3p);
    
    return;
}

void    analyze_read(struct read_sam * p)
{
    int pos, length, flag;
    
    pos = p->Pos;
    length = p->Sequence_length;
    flag = p->Flag;
    
    reads_analyzed++;
    
    if (flag == 0)      // mapped to the + strand, the pos the start of 5' end
    {
        
        if (reference[pos-6] == 'C' && reference[pos+3] == 'G') {flag_5 = 2; S5++;}  // pos-6 is 5' upstream and pos+3 is 5' downstream
        else if (reference[pos-6] == 'C' || reference[pos+3] == 'G') {flag_5 = 1; W5++;}
        else {flag_5 = 0; N5++;}
        
        if (reference[pos+3] == 'G' && reference[pos-6] != 'C') G5++;
        if (reference[pos-6] == 'C' && reference[pos+3] != 'G') C5++;
        
        if (reference[pos-6+length] == 'C' && reference[pos+3+length] == 'G') {flag_3 = 2; S3++;}  // pos-6+length is 3' upstream and pos+3+length is 3' downstream
        else if (reference[pos-6+length] == 'C' || reference[pos+3+length] == 'G') {flag_3 = 1; W3++;}
        else {flag_3 = 0; N3++;}
        
        if (reference[pos+3+length] == 'G' && reference[pos-6+length] != 'C') G3++;
        if (reference[pos-6+length] == 'C' && reference[pos+3+length] != 'G') C3++;
        
        
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
        
        if (reference[pos+5] == 'G' && reference[pos-4] == 'C') {flag_3 = 2; S3++;}  // pos+5 is 3' upstream and pos-4 is 3' downstream
        else if (reference[pos+5] == 'G' || reference[pos-4] == 'C') {flag_3 = 1; W3++;}
        else {flag_3 = 0; N3++;}
        
        if (reference[pos+5] == 'G' && reference[pos-4] != 'C') C3++;
        if (reference[pos-4] == 'C' && reference[pos+5] != 'G') G3++;
        
        
        if (reference[pos+5+length] == 'G' && reference[pos-4+length] == 'C') {flag_5 = 2; S5++;}  // pos+5+length is 5' upstream and pos-4+length is 5' downstream
        else if (reference[pos+5+length] == 'G' || reference[pos-4+length] == 'C') {flag_5 = 1; W5++;}
        else {flag_5 = 0; N5++;}
        
        if (reference[pos-4+length] == 'C' && reference[pos+5+length] != 'G') G5++;
        if (reference[pos+5+length] == 'G' && reference[pos-4+length] != 'C') C5++;
        
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
    
    
    return;
}
    

void rev_com(char* s)
{
    int i, length;
    char d[RADIUS*2+1];
    
    length = RADIUS*2+1;
    
    for (i = 0; i < length - 1; i++)
    {
        if (s[length - 2 - i ] == 'A' ) { d[i] = 'T'; continue;}
        if (s[length - 2 - i ] == 'T' ) { d[i] = 'A'; continue;}
        if (s[length - 2 - i ] == 'G' ) { d[i] = 'C'; continue;}
        if (s[length - 2 - i ] == 'C' ) { d[i] = 'G'; continue;}
    }
    
    d[RADIUS*2] = '\0';
    strcpy (s, d);
}

    


















