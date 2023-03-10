
// 2019 April 2nd

// ***************************************************************************************************************************************
// *    Long_Tail_Analysis_v1.exe input.tsv name_miRNA min_tail_length max_tail_length mode                                              *
// ***************************************************************************************************************************************


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAX_miRNA_SIZE      25              // max length of consensus sequence
#define MAX_tail_SIZE       100              // max length of tail sequence
#define MAX_char_per_line   256             // max size of any line in file


int                 min_tail_length;
int                 max_tail_length;
FILE                *fp_input;
char                name_miRNA[20];             // max size of char = 20!
int                 mode;                   // -1: all; 0: tail-only; other positive number: trimed
int                 tail_composition_table[5];  // 0:total events; 1:A, 2:T 3: G, 4 C;
int                 tail_only_table[5];
int                 tail_trim_table[5];

void    analyze_tail_composition (int*, int, char*, int);

int main(int argc, char* argv[])		// three arguments allowed
{
	
    char NAME[20];
    int READS;
    int LEN_TRIM;
    int LEN_TAIL;
    char SEQ_TAIL[MAX_tail_SIZE];
    
    char                ch;
    char                line[MAX_char_per_line];
    
    int                 num_total_mirna_reads=0;
    int                 num_total_tailed=0;
    int                 num_total_tailed_analyzed=0;
    int                 num_total_tailed_only=0;
    int                 num_total_tailed_trim=0;
    
    
    if (argc != 6)
	{
		fprintf(stderr, "\n Long_Tail_Analysis_v1.exe input.tsv name_miRNA min_tail_length max_tail_length mode \n");
		exit (1);
	}
    
    if ( (fp_input = fopen(argv[1], "r")) == NULL )
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[1]);
        exit (1);
    }
    
    strcpy(name_miRNA, argv[2]);
    
    min_tail_length = atoi (argv[3]);
    max_tail_length = atoi (argv[4]);
    
    mode = atoi (argv[5]);
    
    
    while  ((ch = fgetc(fp_input)) != '\n')                             // skip the first line
    {
        ;
    }
    
    while(fgets(line, MAX_char_per_line, fp_input) != NULL)
    {
       
        
        if (line[0] == 'M') continue;                                      // skip the title line
        sscanf(line, "%s%*s%*d%d%*s%d%d%s", NAME, &READS, &LEN_TRIM, &LEN_TAIL, SEQ_TAIL);
        
 
        if(strcmp(name_miRNA, NAME) != 0) continue;
        num_total_mirna_reads += READS;
        
        if (LEN_TAIL == 0) continue;
        num_total_tailed += READS;
        
        if (LEN_TAIL < min_tail_length || LEN_TAIL > max_tail_length) continue;
        num_total_tailed_analyzed += READS;
        
        analyze_tail_composition (tail_composition_table, LEN_TAIL, SEQ_TAIL, READS);
        if (LEN_TRIM >0)
        {
            analyze_tail_composition (tail_trim_table, LEN_TAIL, SEQ_TAIL, READS);
            num_total_tailed_trim += READS;
        }
        else
        {
            analyze_tail_composition (tail_only_table, LEN_TAIL, SEQ_TAIL, READS);
            num_total_tailed_only += READS;
        }
    }
        
        // output starting here...
        
     
        printf("\n");
        printf("\n");
        printf("%s\t%d\n", name_miRNA, num_total_mirna_reads);
        printf("num_total_tailed\t%d\n",num_total_tailed );
        printf("num_total_tailed between %d and %d is %d\n",min_tail_length, max_tail_length, num_total_tailed_analyzed);
    
        printf("\n");
        printf("\n");
        printf("Mode\tA\tU\tG\tC\tAve_Tail_len\tReads\n");
        printf("All reads\t%.3f\t%.3f\t%.3f\t%.3f\t%.2f\t%d\n", (float)tail_composition_table[1]/(float)tail_composition_table[0], (float)tail_composition_table[2]/(float)tail_composition_table[0], (float)tail_composition_table[3]/(float)tail_composition_table[0], (float)tail_composition_table[4]/(float)tail_composition_table[0], (float)tail_composition_table[0]/(float)num_total_tailed_analyzed, num_total_tailed_analyzed);
        printf("Tail_only\t%.3f\t%.3f\t%.3f\t%.3f\t%.2f\t%d\n", (float)tail_only_table[1]/(float)tail_only_table[0], (float)tail_only_table[2]/(float)tail_only_table[0], (float)tail_only_table[3]/(float)tail_only_table[0], (float)tail_only_table[4]/(float)tail_only_table[0], (float)tail_only_table[0]/(float)num_total_tailed_only, num_total_tailed_only);
        printf("Trim_Tail\t%.3f\t%.3f\t%.3f\t%.3f\t%.2f\t%d\n", (float)tail_trim_table[1]/(float)tail_trim_table[0], (float)tail_trim_table[2]/(float)tail_trim_table[0], (float)tail_trim_table[3]/(float)tail_trim_table[0], (float)tail_trim_table[4]/(float)tail_trim_table[0], (float)tail_trim_table[0]/(float)num_total_tailed_trim, num_total_tailed_trim);
        
       
        
        fclose(fp_input);
        
        return 0;
        
    }


void    analyze_tail_composition (int* table, int LEN_TAIL, char* SEQ_TAIL, int READS)
{
    int    i;
    
    
    for(i=0; i < LEN_TAIL; i++)
    {
        table[0] += READS;
        if (SEQ_TAIL[i] == 'A' )        table[1] += READS;
        else if (SEQ_TAIL[i] == 'T' )   table[2] += READS;
        else if (SEQ_TAIL[i] == 'G' )   table[3] += READS;
        else if (SEQ_TAIL[i] == 'C' )   table[4] += READS;
        else    printf("something is wrong!\n");
    }
}


    /*
    fprintf(fp_output,"Relative Position to the End\tTrimming-only\tTailed\tMono_tailing\tOligo_tailing\tPercentage of Mono-tailing\tPercentage of oligo-tailing\n");
    for(i = upper_trim_size; i >= 0; i--)
    {
        fprintf(fp_output,"%d\t", 0 - i);
        fprintf(fp_output,"%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\t",(float) info[i][1] / (float)num_total_reads*100, (float)(info[i][2] + info[i][3])/(float)num_total_reads*100, (float) info[i][2]/(float)num_total_reads*100, (float) info[i][3]/(float)num_total_reads*100);
        fprintf(fp_output,"%.2f%%\t%.2f%%\n",(float)info[i][2]/(float)(info[i][2] + info[i][3])*100, (float)info[i][3]/(float)(info[i][2] + info[i][3])*100);
    }
    
    fprintf(fp_output,"\n");
    fprintf(fp_output,"\n");
    
    fprintf(fp_output,"Relative Position to the End\tTail-A\tTail-T\tTail-G\tTail-C\tTotal_Mono_tailing_event\t\tRelative Position to the End\tTail-A\tTail-T\tTail-G\tTail-C\tAve_Oligo_Tail_Length\tTotal_Oligo_tailing_event\n");
        
        
        
        
        
        
        
        
        
        
        
        
        
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
            
            if (SEQ_TAIL[1] == 'A' )        A[LEN_TRIM][0] += READS;        // save mono-tail info in the 0 position colume of the tail_composition_table
            else if (SEQ_TAIL[1] == 'T' )   T[LEN_TRIM][0] += READS;
            else if (SEQ_TAIL[1] == 'G' )   G[LEN_TRIM][0] += READS;
            else if (SEQ_TAIL[1] == 'C' )   C[LEN_TRIM][0] += READS;
            else    printf("something is wrong!\n");
        }
        else                                                        // oligo_tailing event
        {
            info[LEN_TRIM][3] += READS;
            info[LEN_TRIM][4] += (READS*LEN_TAIL);
            
            for(i=1; SEQ_TAIL[i] != '\"'; i++)
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
    
    for (i=1; i<= upper_trim_size; i++)
    {
        num_total_trimmed_only += info[i][1];
        num_total_trimmed_and_tailed += (info[i][2] + info[i][3]) ; //info[i][2] - mono tailing, info[i][3] - oligo-tailing
    }
    
    num_total_trimmed = num_total_trimmed_only + num_total_trimmed_and_tailed;
    num_total_tailed  = num_total_tailed_only + num_total_trimmed_and_tailed;
    
    for (i=0; i<= upper_trim_size; i++)
    {
        num_total_mono_tailed   +=  info[i][2];
        num_total_oligo_tailed  +=  info[i][3];
        total_oligo_tail_events +=  info[i][4];
    }
    
    
    
    for(i = upper_trim_size; i >= 0; i--)
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
    
    for (i=0; i<=upper_trim_size; i++)
    {
        fprintf(fp_output,"\t\t\tPosition %d\t\t\t", 0 - i);
    }
    
    fprintf(fp_output,"\t\t\tAll\t\t\t");
    fprintf(fp_output,"\n");
    
    for (i=0; i<=upper_trim_size+1; i++)
    {
        fprintf(fp_output,"\tA\tT\tG\tC\tTotal\t");
    }
    fprintf(fp_output,"\n");
    
    
    for(j = 1; j<= upper_tail_size; j++)                // j point to the position of tail
    {
        
        total_A=0;
        total_T=0;
        total_G=0;
        total_C=0;
        total_pos=0;
        
        for(i=0; i<=upper_trim_size; i++)
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
            
    


long int sum(long int* q )            // p points to N[pos][0]
{
    long int    total=0;
    
    long int*       p;
    
    for (p = q+1; p <= q + upper_tail_size; p++)
    {
        total += *p;
    }
    
    return total;
}
*/
    
