//
//  str_search_mis.c
//  
//
//  Created by Shuo Gu on 11/13/14.
//
//  2014-Dec-18 v2 : Nothing really changed. just updated the strcom_mis to v2. not working now.

#include <stdio.h>
#include <stdlib.h>

char*   strstr_mis (char* s1, char* s2, int n );     /* Returns a pointer to the first occurrence of str2 in str1, or a null pointer if str2 is not part of str1. */
                                                     /* n is the mismatch allowrance */
                                                     /* strstr_mis search the whole s1, and report the FIRST matching of s2 with the LOWEST mismatch number */


int strcmp_mis ( char* s1, char* s2, int n);

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

int strcmp_mis ( char* s1, char* s2, int n)     /* s1 and s2 are the pointers to the strings to compare, n is the numbr of mismatches allowed */
                                                /* for return, -1 - no match, mismatches are more than n,  m>=0 - remaining mismatches allowrance */
                                                /* s1 and s2 both need a \0 at the the ends */
{
    while ( n >= 0 )
    {
        if ( *s1 - *s2 )    n--;
        if ( *s1 )          s1++;
        if ( *s2 )          s2++;
        if ((*s1 == '\0') && (*s2 == '\0')) return n;
    }
    return n;
}
char*   strstr_mis (char* s1, char* s2, int n )     /* Returns a pointer to the first occurrence of str2 in str1, or a null pointer if str2 is not part of str1. */
                                                    /* n is the mismatch allowrance */
                                                    /* strstr_mis search the whole s1, and report the FIRST matching of s2 within the lowest mismatch number */

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
