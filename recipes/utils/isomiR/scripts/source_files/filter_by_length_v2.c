//v2 2015 April 2nd.
// changed the file input and output -- from files


// input is fastq file, output is reads within certain size range
 

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define Bool		int
#define	TRUE		1
#define	FAIL        0
#define PASS        1

#define	LEN_NAME	175
#define	LEN_SEQ		175
#define	LEN_TAG		120


// date structure is here........

struct fastq
{
	char        name[LEN_NAME+1];
	char        sequence[LEN_SEQ+1];
	char		tag[LEN_TAG+1];
	char		score[LEN_SEQ+1];
};

// Function prototypes are here......

struct fastq*	read_fastq(struct fastq*);

void	print_fastq(const struct fastq*);

void	print_fasta( const struct fastq*);


// Global varables are here...

int 	num_read = 0;
float     num_pass = 0;

FILE *fp_input = NULL;
FILE *fp_output = NULL;

// main start from here....

int main(int argc, char* argv[])		// argv[1] and argv[2] will be the min and max of length of reads interested argv[3] will be input file and argv[4] will be output file
{
	if (argc != 5)
	{
		fprintf(stderr, "\n filter_by_length.exe min_length max_length input.fastq output.fastq \n\n");
		return 0;
	}
    
    int min, max, length;
    
    min = atoi(argv[1]);
    max = atoi(argv[2]);
    
    if ( (fp_input = fopen(argv[3], "r")) == NULL )                     // open input file, making sure it is sucessful..
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[3]);
        exit (1);
    }
    fp_output = fopen(argv[4], "w");                      // creat output file.. if file already exits, clear and replace the existing one...

    
    if (max < min)
    {
        fprintf(stderr, "\n mininal length should be shorter or equal to maximal length\n");
        return 0;
    }
    
    struct fastq	read_temp;
	struct fastq 	*p_read = &read_temp;		/* p_read point to the current read */
    
	while (read_fastq(p_read) != NULL )
	{
		num_read++;
        length = strlen(p_read->sequence);
        if( length >= min && length <= max) {num_pass++; print_fastq(p_read);}
	}
	
	fclose(fp_input);
    fclose(fp_output);
    
	fprintf(stdout, "\n");
	fprintf(stdout, "Total number of reads:	%15d\nNumber of reads passing filter ( between %d and %d): %15.0f, that is %.2f%% of total reads\n\n", num_read, min, max, num_pass, num_pass*100/num_read);

	return 0;	
}

struct fastq* read_fastq(struct fastq* read)
{
	char *pos;
    
    if (fgets(read->name, LEN_NAME, fp_input) == NULL) return NULL;
    if ((pos=strchr(read->name, '\n')) != NULL)           *pos = '\0';
	
    fgets(read->sequence, LEN_SEQ, fp_input);
	if ((pos=strchr(read->sequence, '\n')) != NULL)           *pos = '\0';
    
    fgets(read->tag, LEN_TAG, fp_input);
    if ((pos=strchr(read->tag, '\n')) != NULL)           *pos = '\0';
    
	fgets(read->score, LEN_SEQ, fp_input);
    if ((pos=strchr(read->score, '\n')) != NULL)           *pos = '\0';
	
	return read;
}

void print_fastq(const struct fastq* read)
{
	fprintf(fp_output, "%s", read->name);
    fputc( '\n', fp_output);
	fprintf(fp_output, "%s", read->sequence);
    fputc( '\n', fp_output);
	fprintf(fp_output, "%s", read->tag);
    fputc( '\n', fp_output);
	fprintf(fp_output, "%s", read->score);
    fputc( '\n', fp_output);
}

void print_fasta(const struct fastq* read)
{
	putchar('>');
	printf("%s\n", read->name);
	printf("%s\n", read->sequence);
}
