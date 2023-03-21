/*     Original version written by:
 * 	Kevin J. Bowers, Ph.D.
 */

#ifndef MODULE_COMMON
#define MODULE_COMMON

#include <stddef.h> /* For size_t and NULL */

#define BEGIN_PROTOTYPES
#define END_PROTOTYPES

#define BEGIN_PRIMITIVE do
#define END_PRIMITIVE   while(0)

/* The below is an evil preprocessor trick mostly useful for converting the
 * compiler given __LINE__ macro from a number into a string
 */
#define STRINGIFY(x)#x
#define EXPAND_AND_STRINGIFY(x)STRINGIFY(x)

#define CONCAT3(a,b,c)a/**/b/**/c


/*************************************
 * Utilities for descriptive logging *
 *************************************/

/* The argument list to the MESSAGE, WARNING and ERROR macros need to be
 * enclosed in double parens. Example: Assume at file.c line 12 is:
 *
 *   WARNING(("This is a int, %i, and this is a double, %g",1,1.));
 *
 * The following would be printed to the logging facility:
 *   file.c(12): WARNING
 *           This is an int, 1, and this is a double, 1.
 */
#define CHECKPOINT() \
  log_printf(__FILE__"("EXPAND_AND_STRINGIFY(__LINE__)"): Checkpoint\n")

#define MESSAGE(a) BEGIN_PRIMITIVE {                          \
  log_printf(__FILE__"("EXPAND_AND_STRINGIFY(__LINE__)"): "); \
  log_printf a;                                               \
  log_printf("\n");                                           \
} END_PRIMITIVE
    
#define WARNING(a) BEGIN_PRIMITIVE {                                     \
  log_printf(__FILE__"("EXPAND_AND_STRINGIFY(__LINE__)"): WARNING\n\t"); \
  log_printf a;                                                          \
  log_printf("\n");                                                      \
} END_PRIMITIVE

#define ERROR(a) BEGIN_PRIMITIVE {                                       \
  log_printf(__FILE__"("EXPAND_AND_STRINGIFY(__LINE__)"): ERROR\n\t");   \
  log_printf a;                                                          \
  log_printf("\n");                                                      \
} END_PRIMITIVE


/********************************************
 * Utilities for descriptive error messages *
 ********************************************/

typedef const char *error_msg;

#define SUCCESS        ((error_msg)NULL)
#define ERROR_MSG(str) ((error_msg)__FILE__"("EXPAND_AND_STRINGIFY(__LINE__)"): "str)

/* For functions that return error messages, if the function fails, print
   the error message and exit the program */

#define TRY(a) BEGIN_PRIMITIVE { \
  error_msg try_error = (a);     \
  if( try_error ) {              \
    ERROR((try_error));          \
    common_exit();               \
  }                              \
} END_PRIMITIVE


/*******************************************
 * Utilities for dynamic memory allocation *
 *******************************************/

#define MALLOC( p, n ) common_realloc( (void *)&(p), NULL, sizeof(*(p))*(n) )

#define REALLOC( p, mem, n ) \
  common_realloc( (void *)&(p), (void *)(mem), sizeof(*(p))*(n) )

#define FREE( mem ) common_realloc( (void *)&(mem), (void *)(mem), 0 )


#define CLEAR(ptr,n)     memset( (ptr), 0, sizeof(*(ptr))*(n) )

#define COPY(dest,src,n) do {                                           \
    if ((dest)!=(src)) memcpy( (dest), (src), sizeof(*(src))*(n) );     \
  } while(0)


/***********************
 * Function prototypes *
 ***********************/

BEGIN_PROTOTYPES

extern void
log_printf( const char *fmt, ... );

extern void
common_exit( void );

extern void
common_realloc( void * p_ref, /* Use &p */
                void * mem,   /* NULL or a return from previous call */
                size_t size );

extern double
wallclock(void);

END_PROTOTYPES

#endif /* MODULE_COMMON */
