/* read ONE fastq read from stdin and store into a structure, then print into stdout*/


#include <stdio.h>
#include <stdlib.h>

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
void	print_fastq(const struct fastq*);
void print_seq(const struct fastq* );

int main(void)
{
	struct fastq	read_temp;

	for(;;)
	{	

		if (read_fastq(&read_temp) == NULL) return 0;
		
		print_fastq (&read_temp);	

	}
	
}

char* read_fastq(struct fastq* read)
{
	if (gets(read->name) == NULL) return NULL;
	gets(read->sequence);
	gets(read->third);
	gets(read->score);
	
	return read->name;
}

void print_fastq(const struct fastq* read)
{
	printf("%s\n", read->name);
	printf("%s\n", read->sequence);
	printf("%s\n", read->third);
	printf("%s\n", read->score);
}

void print_seq(const struct fastq* read)
{
	printf("%s\n", read->sequence);
}



	