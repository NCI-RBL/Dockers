
// July-28-2018
// take input file txt (line of reads), the cleavage site sequences
// calculate the percentage of -6C and +3G

 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define	LEN_SEQ		175

int main(int argc, char* argv[])		/* three arguments: 1)min 2)input 3) output */
{
	
    FILE                    *fp_input;
	int                     num_total_read = 0;
    int                     both = 0, only_G = 0, only_C = 0, neither = 0;
    char                    read[LEN_SEQ];             // hold each read
    
 //program start here...
    
    if (argc != 2)
	{
		fprintf(stderr, "\n analyze_lan_cleavage_site.exe input.txt \n");
		return 0;
	}
    
    if ( (fp_input = fopen(argv[1], "r")) == NULL )                     // open input file, making sure it is sucessful..
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[1]);
        exit (1);
    }
    
    while (fgets(read, LEN_SEQ, fp_input) != NULL)
    {
        num_total_read++;
        
        if (read[4] == 'C' && read[13] == 'G') both++;
        else if (read[13] == 'G') only_G++;
        else if (read[4] == 'C' ) only_C++;
        else neither++;
    }
    
    fprintf(stdout, "\nTotal number of reads:	%d\n", num_total_read);
    fprintf(stdout, "\n");
    fprintf(stdout, "both +3G and -6C:\t %d\t%d%%\n", both, both*100/num_total_read);
    fprintf(stdout, "either +3G or -6C:\t %d\t%d%%\n", num_total_read - neither, (num_total_read - neither)*100/num_total_read);
    fprintf(stdout, "neither          \t %d\t%d%%\n", neither, neither*100/num_total_read);
    fprintf(stdout, "\n");
    fprintf(stdout, "only +3G         \t %d\t%d%%\n", only_G, only_G*100/num_total_read);
    fprintf(stdout, "only -6C         \t %d\t%d%%\n", only_C, only_C*100/num_total_read);
    
	
    fclose (fp_input);
	return 0;	
}




