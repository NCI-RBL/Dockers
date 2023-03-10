// 2017 Aug 2nd

// ***************************************************************************************************************************************
// *    Tail_Trim_analysis-v1.exe input.tsv output.txt                                                    *
// ***************************************************************************************************************************************

// take Kevin's analysis result (one miRNA per time) as input, perform additional analysis (position-based-information)
// Start from line with '>', put all lines to output file while measure the length of consensus sequence.
// replace Sequence part to new analysis results
// only one miRNA in v1.

// 2017 Aug 8th v2
// separate mono_tail with oligo_tails.
// changed the max_trim to 7 and max_tail to 10.
// skip old analysis results, instead, put all reads info into the new analysis file.
// change the way we analyze the tail composition - all tails aligne at the beginning!
// add length distribution!

// 2017 Aug 11 v3
// make UPPER_trim_SIZE 10 and UPPER_tail_SIZE 11 for the sake of all PAZ data
// add A/T/G/C percentage of 1)mono-tail 2)oligo-tail 3)overall at the end of analysis file


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAX_miRNA_SIZE      25              // max length of consensus sequence
#define MAX_tail_SIZE       40              // max length of tail sequence
#define MAX_char_per_line   256             // max size of any line in file
#define UPPER_trim_SIZE      10              // max length of consensus sequence
#define UPPER_tail_SIZE      11              // max length of tail sequence
#define LENGTH_DISTRIBUTION_UPPER_SIZE      27              
#define LENGTH_DISTRIBUTION_LOWER_SIZE      15



long int sum(long int* p );            // p points to N[pos][0]

int main(int argc, char* argv[])		// three arguments allowed
{
	
    FILE                *fp_input, *fp_output;
    
    char                ch;
    char                line[MAX_char_per_line];
    char                consensus[MAX_miRNA_SIZE];
    int                 size_of_consensus;
    int                 i,j;
    
    long int            info[UPPER_trim_SIZE +1][5] = {0};        // 0-total, 1- No tail, 2- mono tail, 3 - oligo-tail , 4 - number of oligo-event (#3*tail-length)
    long int            A[UPPER_trim_SIZE +1][UPPER_tail_SIZE +1] = {0};        //[pos,0] save mono-event!
    long int            T[UPPER_trim_SIZE +1][UPPER_tail_SIZE +1] = {0};
    long int            G[UPPER_trim_SIZE +1][UPPER_tail_SIZE +1] = {0};
    long int            C[UPPER_trim_SIZE +1][UPPER_tail_SIZE +1] = {0};
    long int            size_distribution[MAX_miRNA_SIZE + MAX_tail_SIZE +1] = {0};
    
    int                 LEN_READ, READS, LEN_TRIM, LEN_TAIL, VAR_5P;
    char                SEQ_TAIL[MAX_tail_SIZE +1];
    
    int                 num_total_reads=0;
    int                 num_skipped_reads=0;
    int                 num_total_isomir=0;
    int                 num_total_trimmed_only=0;
    int                 num_total_trimmed_and_tailed=0;
    int                 num_total_tailed_only=0;
    int                 num_total_trimmed=0;
    int                 num_total_tailed=0;
    int                 num_total_mono_tailed=0;
    int                 num_total_oligo_tailed=0;
    long int            total_oligo_tail_events=0;
    
    long int            total_A=0;
    long int            total_T=0;
    long int            total_G=0;
    long int            total_C=0;
    long int            total_pos=0;
    long int            total=0;
    long int            total_low=0;
    long int            total_high=0;
    
    long int            total_mono_A = 0;
    long int            total_mono_T = 0;
    long int            total_mono_G = 0;
    long int            total_mono_C = 0;
    long int            total_oligo_A = 0;
    long int            total_oligo_T = 0;
    long int            total_oligo_G = 0;
    long int            total_oligo_C = 0;
    long int            total_mono = 0;
    long int            total_oligo = 0;
    

    
    
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
    
    
    while(fgets(line, MAX_char_per_line, fp_input) != NULL)
    {
        if (line[0] == '[') break;                                      // reach the end of the reads section -> stop
        sscanf(line, "%*s%d%d%*s%d%d%s%d", &LEN_READ, &READS, &LEN_TRIM, &LEN_TAIL, SEQ_TAIL, &VAR_5P);
        
        if(VAR_5P != 0) { num_skipped_reads += READS; continue; }                                                                       // skip anything that do not start as where it should be..
        if(LEN_TRIM > UPPER_trim_SIZE || LEN_TAIL > UPPER_tail_SIZE) { num_skipped_reads += READS; continue; }                          // skip anything that has too long tail or trim.
        if(LEN_READ > MAX_miRNA_SIZE + MAX_tail_SIZE) { num_skipped_reads += READS; continue; }                                         // skip anything that is too long.
        
        num_total_reads += READS;
        size_distribution[LEN_READ] += READS;
        fputs(line, fp_output);
        
        info[LEN_TRIM][0] += READS;
        if (LEN_TAIL == 0)                                          // trimming only - No tail
        {
            info[LEN_TRIM][1] += READS;
            continue;
        }
        else if (LEN_TAIL == 1)                                     // momo-tailing event - tail 1
        {
            info[LEN_TRIM][2] += READS;
            
            if (SEQ_TAIL[0] == 'A' )        A[LEN_TRIM][0] += READS;        // save mono-tail info in the 0 position colume of the tail_composition_table
            else if (SEQ_TAIL[0] == 'T' )   T[LEN_TRIM][0] += READS;
            else if (SEQ_TAIL[0] == 'G' )   G[LEN_TRIM][0] += READS;
            else if (SEQ_TAIL[0] == 'C' )   C[LEN_TRIM][0] += READS;
            else    printf("something is wrong!\n");
        }
        else                                                        // oligo_tailing event
        {
            info[LEN_TRIM][3] += READS;
            info[LEN_TRIM][4] += (READS*LEN_TAIL);
            
            for(i=0; SEQ_TAIL[i] != '\0'; i++)
            {
                if (SEQ_TAIL[i] == 'A' )        A[LEN_TRIM][i+1] += READS;
                else if (SEQ_TAIL[i] == 'T' )   T[LEN_TRIM][i+1] += READS;
                else if (SEQ_TAIL[i] == 'G' )   G[LEN_TRIM][i+1] += READS;
                else if (SEQ_TAIL[i] == 'C' )   C[LEN_TRIM][i+1] += READS;
                else    printf("something is wrong!\n");
            }
            
            if (i < 2) printf("something is wrong!\n");            // make sure it is indeed oligo-tailing
        }
    
        
    }
    
    // calculate some statistical info...
    
    
    num_total_isomir = num_total_reads - info[0][1];        //info[0][1] is no trim and no tail -> reads without modification..
    num_total_tailed_only = info[0][2] + info[0][3];
    
    for (i=1; i<= UPPER_trim_SIZE; i++)
    {
        num_total_trimmed_only += info[i][1];
        num_total_trimmed_and_tailed += (info[i][2] + info[i][3]) ; //info[i][2] - mono tailing, info[i][3] - oligo-tailing
    }
    
    num_total_trimmed = num_total_trimmed_only + num_total_trimmed_and_tailed;
    num_total_tailed  = num_total_tailed_only + num_total_trimmed_and_tailed;
    
    for (i=0; i<= UPPER_trim_SIZE; i++)
    {
        num_total_mono_tailed   +=  info[i][2];
        num_total_oligo_tailed  +=  info[i][3];
        total_oligo_tail_events +=  info[i][4];
    }
    
    // output starting here...
    
    fprintf(fp_output,"\n");
    fprintf(fp_output,"\n");
    fprintf(fp_output,"\n");
    fprintf(fp_output,"new analysis start here.. \n");
    fprintf(fp_output,"reads skipped in this analysis\t%d\n",num_skipped_reads);
    fprintf(fp_output,"total-reads\t%d\n",num_total_reads);
    fprintf(fp_output,"num_total_isomir\t%d\t%.2f%%\n",num_total_isomir,(float)num_total_isomir / (float)num_total_reads*100 );
    fprintf(fp_output,"trimming-only\t%.2f%%\n",(float)num_total_trimmed_only / (float)num_total_reads*100);
    fprintf(fp_output,"tailing-only\t%.2f%%\n",(float)num_total_tailed_only / (float)num_total_reads*100);
    fprintf(fp_output,"trimming-and-talling\t%.2f%%\n",(float)num_total_trimmed_and_tailed / (float)num_total_reads*100);
    fprintf(fp_output,"trimming-total\t%.2f%%\n",(float)num_total_trimmed / (float)num_total_reads*100);
    fprintf(fp_output,"tailing-total\t%.2f%%\n",(float)num_total_tailed / (float)num_total_reads*100);
    fprintf(fp_output,"\n");
    fprintf(fp_output,"Mono-tailing reads\t%.2f%%\t%.2f%%\n",(float)num_total_mono_tailed / (float)num_total_reads*100, (float)num_total_mono_tailed / (float)num_total_tailed*100);
    fprintf(fp_output,"Oligo-tailing reads\t%.2f%%\t%.2f%%\n",(float)num_total_oligo_tailed / (float)num_total_reads*100, (float)num_total_oligo_tailed / (float)num_total_tailed*100);
    fprintf(fp_output,"Average Oligo-tail Length\t%.2f\n",(float) total_oligo_tail_events / (float)num_total_oligo_tailed);
    
    
    fprintf(fp_output,"Relative Position to the End\tTrimming-only\tTailed\tMono_tailing\tOligo_tailing\tPercentage of Mono-tailing\tPercentage of oligo-tailing\n");
    for(i = UPPER_trim_SIZE; i >= 0; i--)
    {
        fprintf(fp_output,"%d\t", 0 - i);
        fprintf(fp_output,"%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\t",(float) info[i][1] / (float)num_total_reads*100, (float)(info[i][2] + info[i][3])/(float)num_total_reads*100, (float) info[i][2]/(float)num_total_reads*100, (float) info[i][3]/(float)num_total_reads*100);
        fprintf(fp_output,"%.2f%%\t%.2f%%\n",(float)info[i][2]/(float)(info[i][2] + info[i][3])*100, (float)info[i][3]/(float)(info[i][2] + info[i][3])*100);
    }
    
    fprintf(fp_output,"\n");
    fprintf(fp_output,"\n");
    
    fprintf(fp_output,"Relative Position to the End\tTail-A\tTail-T\tTail-G\tTail-C\tTotal_Mono_tailing_event\t\tRelative Position to the End\tTail-A\tTail-T\tTail-G\tTail-C\tAve_Oligo_Tail_Length\tTotal_Oligo_tailing_event\n");
    
    
    for(i = UPPER_trim_SIZE; i >= 0; i--)
    {
        fprintf(fp_output,"%d\t", 0 - i);
        
        fprintf(fp_output,"%.2f%%\t",(float) A[i][0] / (float) info[i][2]*100);
        fprintf(fp_output,"%.2f%%\t",(float) T[i][0] / (float) info[i][2]*100);
        fprintf(fp_output,"%.2f%%\t",(float) G[i][0] / (float) info[i][2]*100);
        fprintf(fp_output,"%.2f%%\t",(float) C[i][0] / (float) info[i][2]*100);
        
        fprintf(fp_output,"%ld\t", info[i][2]);
        
        fprintf(fp_output,"\t%d\t", 0 - i);
        
        fprintf(fp_output,"%.2f%%\t",(float) sum(A[i]) / (float) info[i][4]*100);
        fprintf(fp_output,"%.2f%%\t",(float) sum(T[i]) / (float) info[i][4]*100);
        fprintf(fp_output,"%.2f%%\t",(float) sum(G[i]) / (float) info[i][4]*100);
        fprintf(fp_output,"%.2f%%\t",(float) sum(C[i]) / (float) info[i][4]*100);
        
        fprintf(fp_output,"%.2f\t",(float) info[i][4] / (float) info[i][3]);
        fprintf(fp_output,"%ld", info[i][4]);
        fprintf(fp_output,"\n");
        
        total_mono_A += A[i][0];
        total_mono_T += T[i][0];
        total_mono_G += G[i][0];
        total_mono_C += C[i][0];
        
        total_oligo_A += sum(A[i]);
        total_oligo_T += sum(T[i]);
        total_oligo_G += sum(G[i]);
        total_oligo_C += sum(C[i]);
    }
    
    fprintf(fp_output,"\n");

    //print more info of oligo-tail composition by position
    
    for (i=0; i<=UPPER_trim_SIZE; i++)
    {
        fprintf(fp_output,"\t\t\tPosition %d\t\t\t", 0 - i);
    }
    
    fprintf(fp_output,"\t\t\tAll\t\t\t");
    fprintf(fp_output,"\n");
    
    for (i=0; i<=UPPER_trim_SIZE+1; i++)
    {
        fprintf(fp_output,"\tA\tT\tG\tC\tTotal\t");
    }
    fprintf(fp_output,"\n");
    
    
    for(j = 1; j<= UPPER_tail_SIZE; j++)                // j point to the position of tail
    {
        
        total_A=0;
        total_T=0;
        total_G=0;
        total_C=0;
        total_pos=0;
        
        for(i=0; i<=UPPER_trim_SIZE; i++)
        {
            total = A[i][j] + T[i][j] + G[i][j] + C[i][j];
            fprintf(fp_output,"%d\t", j );
            fprintf(fp_output,"%.2f%%\t",(float) A[i][j] / (float) total*100);
            fprintf(fp_output,"%.2f%%\t",(float) T[i][j] / (float) total*100);
            fprintf(fp_output,"%.2f%%\t",(float) G[i][j] / (float) total*100);
            fprintf(fp_output,"%.2f%%\t",(float) C[i][j] / (float) total*100);
            fprintf(fp_output,"%ld\t", total);          // for tail 1 and tail 2, total should be the same as oligo-tail info[i][3], for tail 3 and later, the number will do down a lot..
            
            total_pos += total;
            total_A += A[i][j];
            total_T += T[i][j];
            total_G += G[i][j];
            total_C += C[i][j];
        }
        
        fprintf(fp_output,"%d\t", j );
        fprintf(fp_output,"%.2f%%\t",(float) total_A / (float) total_pos*100);
        fprintf(fp_output,"%.2f%%\t",(float) total_T / (float) total_pos*100);
        fprintf(fp_output,"%.2f%%\t",(float) total_G / (float) total_pos*100);
        fprintf(fp_output,"%.2f%%\t",(float) total_C / (float) total_pos*100);
        fprintf(fp_output,"%ld\t", total_pos);
        fprintf(fp_output,"\n");
    }
    

    fprintf(fp_output,"\n");
    
    // print length distribution...
    
    fprintf(fp_output,"Length Distribution\n");
    fprintf(fp_output,"Position\t<%d\t", LENGTH_DISTRIBUTION_LOWER_SIZE);
    
    for(i= LENGTH_DISTRIBUTION_LOWER_SIZE; i <= LENGTH_DISTRIBUTION_UPPER_SIZE; i++)
    {
        fprintf(fp_output,"%d\t", i);
    }
    
    fprintf(fp_output,">%d\t", LENGTH_DISTRIBUTION_UPPER_SIZE);
    fprintf(fp_output,"\n");
    
    
    fprintf(fp_output,"Abundance\t");
            
    
    for (i = 1; i < LENGTH_DISTRIBUTION_LOWER_SIZE; i++)
    {
        total_low += size_distribution[i];
    }
    fprintf(fp_output,"%ld\t", total_low);
    
    for(i= LENGTH_DISTRIBUTION_LOWER_SIZE; i <= LENGTH_DISTRIBUTION_UPPER_SIZE; i++)
    {
        fprintf(fp_output,"%ld\t", size_distribution[i]);
    }
    
    for (i = LENGTH_DISTRIBUTION_UPPER_SIZE +1; i< MAX_miRNA_SIZE + MAX_tail_SIZE; i++)
    {
        total_high += size_distribution[i];
    }
    fprintf(fp_output,"%ld\t", total_high);
    fprintf(fp_output,"\n");
    
    
            
    fprintf(fp_output,"Percentage\t");
            
    fprintf(fp_output,"%.2f%%\t", (float)total_low/(float)num_total_reads*100);
                    
    for(i= LENGTH_DISTRIBUTION_LOWER_SIZE; i <= LENGTH_DISTRIBUTION_UPPER_SIZE; i++)
    {
        fprintf(fp_output,"%.2f%%\t", (float)size_distribution[i]/(float)num_total_reads*100);
    }
                    
    fprintf(fp_output,"%.2f%%\t", (float)total_high/(float)num_total_reads*100);
    fprintf(fp_output,"\n");
    fprintf(fp_output,"Tail composition analysis.....\n");
    fprintf(fp_output,"\n");
    fprintf(fp_output,"Mono\tTail-A\tTail-T\tTail-G\tTail-C\tTotal_Mono_tailing_event\t\tOligo\tTail-A\tTail-T\tTail-G\tTail-C\tTotal_Oligo_tailing_event\t\tOverall\tTail-A\tTail-T\tTail-G\tTail-C\tTotal_tailing_event\n");
    total_mono = total_mono_A + total_mono_T + total_mono_G + total_mono_C;
    total_oligo = total_oligo_A + total_oligo_T + total_oligo_G + total_oligo_C;
    fprintf(fp_output,"\t");
    fprintf(fp_output,"%.2f%%\t",(float) total_mono_A / (float) total_mono*100);
    fprintf(fp_output,"%.2f%%\t",(float) total_mono_T / (float) total_mono*100);
    fprintf(fp_output,"%.2f%%\t",(float) total_mono_G / (float) total_mono*100);
    fprintf(fp_output,"%.2f%%\t",(float) total_mono_C / (float) total_mono*100);
    fprintf(fp_output,"%ld\t", total_mono);
            
    fprintf(fp_output,"\t");
    fprintf(fp_output,"\t");
    fprintf(fp_output,"%.2f%%\t",(float) total_oligo_A / (float) total_oligo*100);
    fprintf(fp_output,"%.2f%%\t",(float) total_oligo_T / (float) total_oligo*100);
    fprintf(fp_output,"%.2f%%\t",(float) total_oligo_G / (float) total_oligo*100);
    fprintf(fp_output,"%.2f%%\t",(float) total_oligo_C / (float) total_oligo*100);
    fprintf(fp_output,"%ld\t", total_oligo);
            
    fprintf(fp_output,"\t");
    fprintf(fp_output,"\t");
    fprintf(fp_output,"%.2f%%\t",(float) (total_mono_A + total_oligo_A) / (float) (total_mono + total_oligo)*100);
    fprintf(fp_output,"%.2f%%\t",(float) (total_mono_T + total_oligo_T) / (float) (total_mono + total_oligo)*100);
    fprintf(fp_output,"%.2f%%\t",(float) (total_mono_G + total_oligo_G) / (float) (total_mono + total_oligo)*100);
    fprintf(fp_output,"%.2f%%\t",(float) (total_mono_C + total_oligo_C) / (float) (total_mono + total_oligo)*100);
    fprintf(fp_output,"%ld\t", total_mono + total_oligo);
            
    fclose(fp_input);
    fclose(fp_output);
    
    return 0;
    
}


long int sum(long int* q )            // p points to N[pos][0]
{
    long int    total=0;
    
    long int*       p;
    
    for (p = q+1; p <= q + UPPER_tail_SIZE; p++)
    {
        total += *p;
    }
    
    return total;
}

