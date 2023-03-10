//
//  strcmp_mis.c
//  
//
//  Created by Shuo Gu on 11/12/14.
//
//

#include <stdio.h>
#include <stdlib.h>


int strcmp_mis ( char* s1, char* s2, int n);    /* s1 is the pointer to the subject string, s2 substring, n is the numbr of mismatches allowed */
                                                /* for return, -1 - no match, mismatches are more than n, -2 - s1 is shorter then s2, m>=0 - remaining mismatches allowrance */
                                                /* s2 need a \0 at the the end, no such requirment on s1 */

int main (int argc, char* argv[])
{
    int n;
    
    n = atoi(argv[3]);
    
    n = strcmp_mis(argv[1], argv[2], n);
    
    if ( n >=0) printf("\nMatches, the remaining mismatch allowrance is %d\n", n);
    else if ( n == -1) printf("\n Too many mismatches!\n");
    else if ( n == -2) printf("\n s1 is shorter than s2!\n");
    else printf("Error!, n has a value of %d\n", n);
    
    return 0;
    
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