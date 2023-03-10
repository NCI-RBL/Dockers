// 2017 Aug 2nd

// ***************************************************************************************************************************************
// *    Tail_Trim_analysis-v1.exe input.tsv output.txt                                                    *
// ***************************************************************************************************************************************

// take Kevin's analysis result (one miRNA per time) as input, perform additional analysis (position-based-information)
// Start from line with '>', put all lines to output file while measure the length of consensus sequence.
// replace Sequence part to new analysis results
// only one miRNA in v1.


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAX_miRNA_SIZE      25              // max length of consensus sequence
#define MAX_tail_SIZE       40              // max length of tail sequence
#define MAX_char_per_line   256             // max size of any line in file
#define UPPER_trim_SIZE      8              // max length of consensus sequence
#define UPPER_tail_SIZE      20              // max length of tail sequence


long int sum(long int* p );            // p points to N[pos][0]

int main(int argc, char* argv[])		// three arguments allowed
{
	
    FILE                *fp_input, *fp_output;
    
    char                ch;
    char                line[MAX_char_per_line];
    char                consensus[MAX_miRNA_SIZE];
    int                 size_of_consensus;
    int                 i, j, k;
    
    long int            info[MAX_miRNA_SIZE][5] = {0};        // 0-relative position, 1-Trimmed-only, 2-Trimmed_Tailed, 3-Trimmed(1+2), 4-number of total tailing events
    long int            A[MAX_miRNA_SIZE][MAX_tail_SIZE] = {0};
    long int            T[MAX_miRNA_SIZE][MAX_tail_SIZE] = {0};
    long int            G[MAX_miRNA_SIZE][MAX_tail_SIZE] = {0};
    long int            C[MAX_miRNA_SIZE][MAX_tail_SIZE] = {0};
    
    int                 READS, LEN_TRIM, LEN_TAIL, VAR_5P;
    char                SEQ_TAIL[MAX_tail_SIZE];
    
    int                 num_total_reads=0;
    int                 num_total_isomir=0;
    int                 num_total_trimmed_only=0;
    int                 num_total_trimmed_and_tailed=0;
    int                 num_total_tailed_only=0;
    int                 num_total_trimmed=0;
    int                 num_total_tailed=0;
    
    int                 pos;
    long int            total;
    long int            total_A=0;
    long int            total_T=0;
    long int            total_G=0;
    long int            total_C=0;
    long int            total_pos=0;
    
    if (argc != 3)
	{
		fprintf(stderr, "\n Tail_Trim_analysis-v1.exe input.tsv output.txt  \n");
		exit (1);
	}
    
    if ( (fp_input = fopen(argv[1], "r")) == NULL )
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[1]);
        exit (1);
    }
    
    fp_output = fopen(argv[2], "w");
    
    while  ((ch = fgetc(fp_input)) != '>')                             // skip everything before the sign '>'
    {
        ;
    }
    fputc(ch, fp_output);                                               // print '>' into the output file
    fgets(line, MAX_char_per_line, fp_input); fputs(line, fp_output);   // rest of the first line
    fgets(line, MAX_char_per_line, fp_input); fputs(line, fp_output);   // second line (motif)
    fgets(line, MAX_char_per_line, fp_input); fputs(line, fp_output);   // third line (consensus)
    
    sscanf(line, "%*s%s", consensus);
    size_of_consensus = strlen (consensus);
    
    for(i=1; i<15; i++)                                                 // skip 14 lines from input
    {
        fgets(line, MAX_char_per_line, fp_input);
        //fputs(line, fp_output);
    }
    
    //fgets(line, MAX_char_per_line, fp_input);
    //fgets(line, MAX_char_per_line, fp_input);                           // skip two lines
    
    while(fgets(line, MAX_char_per_line, fp_input) != NULL)
    {
        if (line[0] == '[') break;                                      // reach the end of the reads section -> stop
        sscanf(line, "%*s%*d%d%*s%d%d%s%d", &READS, &LEN_TRIM, &LEN_TAIL, SEQ_TAIL, &VAR_5P);
        
        if(VAR_5P != 0) continue;                                       // skip anything that do not start as where it should be..
        if(LEN_TRIM > UPPER_trim_SIZE || LEN_TAIL > UPPER_tail_SIZE) continue;                         // skip anything that has too long tail or trim.
    
        num_total_reads += READS;
        
        pos=size_of_consensus - LEN_TRIM;
        
        if (LEN_TAIL == 0)
        {
            info[pos][1] += READS;
            continue;
        }
        else
        {
            info[pos][2] += READS;
            info[pos][4] += (READS*LEN_TAIL);
        }
        
        // have tail, so need to fill up the tail composition info table...
        
        for(i=0; SEQ_TAIL[i] != '\0'; i++)
        {
            if (SEQ_TAIL[i] == 'A' ) A[pos][i+1] += READS;
            else if (SEQ_TAIL[i] == 'T' ) T[pos][i+1] += READS;
            else if (SEQ_TAIL[i] == 'G' ) G[pos][i+1] += READS;
            else if (SEQ_TAIL[i] == 'C' ) C[pos][i+1] += READS;
            else    printf("something is wrong!\n");
        }
    }
    
    // file-up the rest of info-table..
    
    for (i=1; i<= size_of_consensus; i++)
    {
        info[i][0] = i - size_of_consensus;
        info[i][3] = info[i][1] + info[i][2];
        
        num_total_isomir += info[i][3];
        num_total_trimmed_only += info[i][1];
        num_total_trimmed_and_tailed += info[i][2];
    }
    
    num_total_isomir -= info[size_of_consensus][1];
    num_total_trimmed_only -= info[size_of_consensus][1];
    num_total_trimmed_and_tailed -= info[size_of_consensus][2];
    num_total_tailed_only = info[size_of_consensus][2];
    num_total_trimmed=num_total_trimmed_only+num_total_trimmed_and_tailed;
    num_total_tailed=num_total_tailed_only+num_total_trimmed_and_tailed;
    
    fprintf(fp_output,"\n");
    fprintf(fp_output,"new analysis start here.. \n");
    fprintf(fp_output,"total-reads\t%d\n",num_total_reads);
    fprintf(fp_output,"num_total_isomir\t%d\t%.1f%%\n",num_total_isomir,(float)num_total_isomir / (float)num_total_reads*100 );
    fprintf(fp_output,"trimming-only\t%.1f%%\n",(float)num_total_trimmed_only / (float)num_total_reads*100);
    fprintf(fp_output,"tailing-only\t%.1f%%\n",(float)num_total_tailed_only / (float)num_total_reads*100);
    fprintf(fp_output,"trimming-and-talling\t%.1f%%\n",(float)num_total_trimmed_and_tailed / (float)num_total_reads*100);
    fprintf(fp_output,"trimming-total\t%.1f%%\n",(float)num_total_trimmed / (float)num_total_reads*100);
    fprintf(fp_output,"tailing-total\t%.1f%%\n",(float)num_total_tailed / (float)num_total_reads*100);
    fprintf(fp_output,"\n");
    
    fprintf(fp_output,"Position\tRelative Position to the End\tTrimming-only\tTrimmed_and_Tailed\tAve_Tail_Length\tTail-A\tTail-T\tTail-G\tTail-C\ttotal tailing events\n");
    for(pos=size_of_consensus - UPPER_trim_SIZE; pos <= size_of_consensus; pos++)
    {
        fprintf(fp_output,"%d\t%ld\t", pos, info[pos][0]);
        fprintf(fp_output,"%.2f%%\t%.2f%%\t",(float) info[pos][1] / (float)num_total_reads*100,(float) info[pos][2]/(float)num_total_reads*100);
        fprintf(fp_output,"%.2f\t",(float) info[pos][4] / (float) info[pos][2]);
        
        fprintf(fp_output,"%.2f%%\t",(float) sum(A[pos]) / (float) info[pos][4]*100);
        fprintf(fp_output,"%.2f%%\t",(float) sum(T[pos]) / (float) info[pos][4]*100);
        fprintf(fp_output,"%.2f%%\t",(float) sum(G[pos]) / (float) info[pos][4]*100);
        fprintf(fp_output,"%.2f%%\t",(float) sum(C[pos]) / (float) info[pos][4]*100);
        
        fprintf(fp_output,"%ld", info[pos][4]);
        fprintf(fp_output,"\n");
    }
    
    //print table of tail info..
    
    for (i=0; i<=UPPER_trim_SIZE; i++)
    {
        fprintf(fp_output,"\t\t\tPosition %d\t\t\t", 0-UPPER_trim_SIZE+i);
    }
    fprintf(fp_output,"\t\t\tAll\t\t\t");
    fprintf(fp_output,"\n");
    
    for (i=0; i<=UPPER_trim_SIZE+1; i++)
    {
        fprintf(fp_output,"\tA\tT\tG\tC\tTotal\t");
    }
    fprintf(fp_output,"\n");
    
    for(i = 0 - UPPER_trim_SIZE + 1; i <=  UPPER_tail_SIZE; i++)                  // i point to the position in tail
    {
        fprintf(fp_output,"%d", i );
        
        total_A=0;
        total_T=0;
        total_G=0;
        total_C=0;
        total_pos=0;
        
        for (j = 0 - UPPER_trim_SIZE; (j < i) && (j <= 0); j++)                  // j point to the position of trim
        {
            pos = size_of_consensus + j;
            total = 0;
            //printf("i = %d, j= %d, pos = %d\n", i, j, pos);
            total = A[pos][i-j]+T[pos][i-j]+G[pos][i-j]+C[pos][i-j];
            fprintf(fp_output,"\t%.1f%%\t",(float) A[pos][i-j] / (float) total*100);
            fprintf(fp_output,"%.1f%%\t",(float) T[pos][i-j] / (float) total*100);
            fprintf(fp_output,"%.1f%%\t",(float) G[pos][i-j] / (float) total*100);
            fprintf(fp_output,"%.1f%%\t",(float) C[pos][i-j] / (float) total*100);
            fprintf(fp_output,"%ld\t", total);
            
            total_A += A[pos][i-j];
            total_T += T[pos][i-j];
            total_G += G[pos][i-j];
            total_C += C[pos][i-j];
            total_pos += total;
            
        }
        
        for (k = 0; k < (0-i) + 1; k++)
        {
            fprintf(fp_output,"\t\t\t\t\t\t");
        }
        
        
        fprintf(fp_output,"\t%.1f%%\t",(float) total_A / (float) total_pos*100);
        fprintf(fp_output,"%.1f%%\t",(float) total_T / (float) total_pos*100);
        fprintf(fp_output,"%.1f%%\t",(float) total_G / (float) total_pos*100);
        fprintf(fp_output,"%.1f%%\t",(float) total_C / (float) total_pos*100);
        fprintf(fp_output,"%ld\t", total_pos);
     
        fprintf(fp_output,"\n");
    }

    fprintf(fp_output,"\n");
    
    
    //printf("the content of consensus is %s\n", consensus);
    //printf("the size of consensus is %d\n", size_of_consensus);
    
    //printf("READS is %d\n", READS);
    //printf("LEN_TRIM is %d\n", LEN_TRIM);
    //printf("LEN_TAIL is %d\n", LEN_TAIL);
    //printf("SEQ_TAIL is %s\n", SEQ_TAIL);
    //printf("VAR_5P is %d\n", VAR_5P);
    
    fclose(fp_input);
    fclose(fp_output);
    
    return 0;
    
}


long int sum(long int* q )            // p points to N[pos][0]
{
    long int    total=0;
    int         i;
    long int*       p;
    
    for (p = q+1; p <= q + UPPER_tail_SIZE; p++)
    {
        total += *p;
    }
    
    return total;
}

