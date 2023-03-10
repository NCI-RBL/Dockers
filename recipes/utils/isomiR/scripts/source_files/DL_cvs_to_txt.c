// 2016 Nov 18th

// ***************************************************************************************************************************************
// *    DL_cvs_to_txt.exe input_1.cvs output.txt                                                    *
// ***************************************************************************************************************************************

// used to convert cvs to txt: 1. change ^M (13) to '\n' 2. change ','(44) to '\t'

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAX_LINE_SIZE   100


int main(int argc, char* argv[])		// three arguments allowed
{
	
    FILE                *fp_input, *fp_output;
    
    char                ch;
    
    
    if (argc != 3)
	{
		fprintf(stderr, "\n DL_cvs_to_txt.exe input_1.cvs output.txt  \n");
		exit (1);
	}
    
    if ( (fp_input = fopen(argv[1], "r")) == NULL )
    {
        fprintf(stderr, "Cannot find file %s!\n", argv[1]);
        exit (1);
    }
    
    fp_output = fopen(argv[2], "w");
    
    while  ((ch = fgetc(fp_input)) != EOF)                             // skp the first line for file1
    {
        if (ch == 13) ch = '\n';
        if (ch == 44) ch = '\t';
        fputc(ch, fp_output);
    }
    
    fclose(fp_input);
    fclose(fp_output);
    
    return 0;
    
}
