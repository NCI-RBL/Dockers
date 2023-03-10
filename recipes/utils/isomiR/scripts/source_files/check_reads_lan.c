
// July-23rd, 2016
// This program read fastq file from the stdin, check the composition at specific position of the reads and give summary to stdout
// Check if the +4 position is G? Check if the the -6 position is C?
// How many fit one condition, how many fit two coditions etc..



#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define	LEN_NAME	70
#define	LEN_SEQ		75
	
struct fastq 
{
	char 	name[LEN_NAME+1];
	char 	sequence[LEN_SEQ+1];
	char	third[2];
	char	score[LEN_SEQ+1];
};

char*	read_fastq(struct fastq*);


int main(void)
{
	struct fastq	read_temp;
    char            read[LEN_SEQ+1];
    int             length;
    int             both=0, G4=0, C6=0, neither=0;
    int             number_read=0;
    
	for(;;)
	{	

        if (read_fastq(&read_temp) == NULL) break;
        
        number_read++;
		strcpy(read, read_temp.sequence);
        length = strlen (read);
        
        if (read[3] == 'G' && read[length - 6] == 'C') both++;
        if (read[3] == 'G' && read[length - 6] != 'C') G4++;
        if (read[3] != 'G' && read[length - 6] == 'C') C6++;
        if (read[3] != 'G' && read[length - 6] != 'C') neither++;
	}
	
    printf("total number of reads is: %d\n", number_read);
    printf("both\tG4\tC6\tneither\n");
    
    printf("%4d%%\t%4d%%\t%4d%%\t%4d%%\n", both*100/number_read, G4*100/number_read, C6*100/number_read, neither*100/number_read);
    
    return 0;
}

char* read_fastq(struct fastq* read)
{
	if (gets(read->name) == NULL) return NULL;
    gets(read->sequence);
	gets(read->third);
	gets(read->score);
	
	return read->name;
}




	