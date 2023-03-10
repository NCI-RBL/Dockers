//
//  strcmp_mis.c
//  
//
//  Created by Shuo Gu on 11/12/14.
//
//  2014-Dec-18 v2 : remove special requirements on both s1 and s2.

#include <stdio.h>
#include <stdlib.h>


int strcmp_mis ( char* s1, char* s2, int n);


int main (int argc, char* argv[])
{
    int n;
    
    n = atoi(argv[3]);
    
    n = strcmp_mis(argv[1], argv[2], n);
    
    if ( n >=0) printf("\nMatches, the remaining mismatch allowrance is %d\n", n);
    else if ( n == -1) printf("\n Too many mismatches!\n");
    else printf("Error!, n has a value of %d\n", n);
    
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