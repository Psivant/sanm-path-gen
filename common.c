/* Written by Kevin J. Bowers */

#include <stdarg.h>   /* For va_list, va_start, va_end */
#include <stdio.h>    /* For stderr and vfprintf */
#include <stdlib.h>   /* For malloc, free, realloc and exit */
#include <string.h>   /* For memset, memcpy */
#include <sys/time.h> /* For gettimeofday */

#include "common.h"   /* For function prototypes and NULL */

void
log_printf( const char *fmt, ... ) {
  va_list ap;
  
  /* Check the input arguments */

  if( fmt==NULL ) return;
  
  /* Print the string to the stream */

  va_start(ap,fmt);
  vfprintf(stderr,fmt,ap);
  va_end(ap);

  /* Force string to print before returning */

  fflush(stderr);
}

void
common_exit( void ) {
  exit(1);
}


void
common_realloc( void * mem_ref, /* Use: &p */
                void * mem_old,
                size_t size ) {

# define mem (*((void **)mem_ref))

  /* Check input arguments */

  if( mem_ref==NULL && size!=0 ) ERROR(( "Request would leak" ));

  if( size==0 ) { /* A free request ... mem_ref may be NULL */

    if( mem_ref!=NULL ) mem = NULL;
    if( mem_old!=NULL ) free( mem_old );

  } else if( mem_old==NULL ) { /* A malloc request ... mem_ref is not NULL */

    mem = malloc( size );
    if( mem==NULL ) ERROR(( "Malloc failed" ));

  } else { /* A realloc request ... mem_ref and mem_old are not NULL */
    
    mem = realloc( mem_old, size );
    if( mem==NULL ) ERROR(( "Realloc failed" ));

  }

# undef mem
}


double
wallclock( void ) {
  struct timeval tv;

  gettimeofday( &tv, NULL );
  return tv.tv_sec + tv.tv_usec*1e-6;
}

