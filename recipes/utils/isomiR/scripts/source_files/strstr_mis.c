//
//  str_search_mis.c
//  
//
//  Created by Shuo Gu on 11/13/14.
//
//

#include <stdio.h>
#include <stdlib.h>

char*   strstr_mis (char* s1, char* s2, int n );     /* Returns a pointer to the first occurrence of str2 in str1, or a null pointer if str2 is not part of str1. */
                                                     /* n is the mismatch allowrance */
                                                     /* strstr_mis search the whole s1, and report the FIRST matching of s2 with the LOWEST mismatch number */


int strcmp_mis ( char* s1, char* s2, int n);        /* s1 is the pointer to the subject string, s2 substring, n is the numbr of mismatches allowed */
                                                    /* for return, -1 - no match, mismatches are more than n, -2 - s1 is shorter then s2, m>=0 - remaining mismatches allowrance */
                                                    /* s2 need a \0 at the the end, no such requirment on s1 */

int main (int argc, char* argv[])
{
    int n;
    char * pos;
    
    n = atoi(argv[3]);
    
    pos = strstr_mis(argv[1], argv[2], n);
    
    if ( pos == NULL) printf("\nNo match Found!\n");
    else printf("Found! %s\n", pos);
    
    return 0;
    
}

int strcmp_mis ( char* s1, char* s2, int n)     /* s1 is the pointer to the subject string, s2 substring, n is the numbr of mismatches allowed */
                                                /* for return, -1 - no match, mismatches are more than n, -2 - s1 is shorter then s2, m>=0 - remaining mismatches allowrance */
                                                /* s2 need a \0 at the the end, no such requirment on s1 */
{
    while ( *s2 )
    {
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
