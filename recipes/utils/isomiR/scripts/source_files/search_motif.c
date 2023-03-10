/* this grogram search specific sequence motif in all reads of a fastq file*/
/* and then print out the statistics */
/* search_motif MOTIF number_of_missmatch_allowed < file.fastq */
/* modification made in 11/4/2014 -- 1. allow mismatch 2. tell position info */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#define	LEN_NAME	75
#define	LEN_SEQ		150
	
struct fastq 
{
	char 	name[LEN_NAME+1];
	char 	sequence[LEN_SEQ+1];
	char	third[2];
	char	score[LEN_SEQ+1];
};

char*	read_fastq(struct fastq*);
void	print_fastq(const struct fastq*);
void 	print_seq(const struct fastq* );
int	search(const struct fastq*, const char*, int );   /* first argument is the pointer to read, second is pointer to searching keyword, third is the number of mismatches allowed  */
char*   strstr_mis (char* s1, char* s2, int n );     /* Returns a pointer to the first occurrence of str2 in str1, or a null pointer if str2 is not part of str1. */
                                                        /* n is the mismatch allowrance */
                                                        /* strstr_mis search the whole s1, and report the FIRST matching of s2 with the LOWEST mismatch number */
int strcmp_mis ( char* s1, char* s2, int n);        /* s1 is the pointer to the subject string, s2 substring, n is the numbr of mismatches allowed */
                                                    /* for return, -1 - no match, mismatches are more than n, -2 - s1 is shorter then s2, m>=0 - remaining mismatches allowrance */
                                                    /* s2 need a \0 at the the end, no such requirment on s1 */


int main(int argc, char* argv[])		/* argv[1] will be the pointer to search motif, argv[2] will be the number of missmatch allowed*/
{
	if (argc != 3)
	{
		printf("search_motif MOTIF number_of_missmatch_allowed < file.fastq\n");
		return 0;
	}

	long int 	num_read = 0, no_match = 0, one_match = 0, two_match = 0, more_match = 0;
    long int    num_one_mismatch = 0, num_two_mismatches = 0, num_more_mismatches =0;
	int		motif_count, n;

	struct fastq	read_temp;
	struct fastq 	*p_read = &read_temp;	/* p_read point to the current read */

	n = atoi(argv[2]);
    if ( strlen (argv[1]) <= n) exit (1);  /* too many mismatches allowed */
    
    while (read_fastq(p_read) != NULL )
	{	
		num_read++;
		motif_count = search(p_read, argv[1], n);
		switch (motif_count)
		{
			case 0:	no_match++; break;
			case 1: one_match++; break;
			case 2: two_match++; break;
			default: more_match++;
		}	
	}
	
	printf("Total reads %ld\n", num_read);
	printf("number of reads have no match:              %15ld that is %ld%% of total reads\n", no_match, no_match*100/num_read);
	printf("number of reads have one match:             %15ld that is %ld%% of total reads\n", one_match, one_match*100/num_read);
	printf("number of reads have two matches:           %15ld that is %ld%% of total reads\n", two_match, two_match*100/num_read);
	printf("number of reads have more than two matches: %15ld that is %ld%% of total reads\n", more_match, more_match*100/num_read);
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

int  search(const struct fastq* read, const char* keyword, int n )
{
	const char*	p = read->sequence;
	int	count = 0;

	while ( (p = strstr_mis(p, keyword, n)) != NULL)
	{
		count++;
		p++;
	}

	return count;
}
int strcmp_mis ( char* s1, char* s2, int n)     /* s1 is the pointer to the subject string, s2 substring, n is the numbr of mismatches allowed */
/* for return, -1 - no match, mismatches are more than n, -2 - s1 is shorter then s2, m>=0 - remaining mismatches allowrance */
/* s2 need a \0 at the the end, no such requirment on s1 */
{
    while ( *s2 )
    {
        if (*s1 == '\0')    return -2;
        if ((*s1++ - *s2++) && !(n--)) break;
    }
    return n;
}

char*   strstr_mis (char* s1, char* s2, int n )     /* Returns a pointer to the first occurrence of str2 in str1, or a null pointer if str2 is not part of str1. */
/* n is the mismatch allowrance */
/* strstr_mis search the whole s1, and report the FIRST matching of s2 with the LOWEST mismatch number */

{
    char * pos = NULL;
    int min_mis = -1, mismatch;
    
    while (*s1)
    {
        mismatch = strcmp_mis(s1, s2, n);
        if (mismatch > min_mis) min_mis = mismatch, pos = s1;
        s1++;
    }
    
    return pos;
}








