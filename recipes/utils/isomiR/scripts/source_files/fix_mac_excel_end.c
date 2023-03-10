// 2017 Aug 9

// ***************************************************************************************************************************************
// *    fix_mac_excel_end.exe input.txt output.txt                                                    *
// ***************************************************************************************************************************************

// used to convert txt to txt: change ^M (13) to '\n'

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
    
    while  ((ch = fgetc(fp_input)) != EOF)
    {
        if (ch == 13) ch = '\n';
        fputc(ch, fp_output);
    }
    
    fclose(fp_input);
    fclose(fp_output);
    
    return 0;
    
}
