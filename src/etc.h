#ifndef _ETC_H_
#define _ETC_H_

static void
nofun()
{
}

#ifndef NOT_R_PACKAGE

#include <R.h>
#include <R_ext/Utils.h>

#define mess(...) Rprintf(__VA_ARGS__)

#define nalloc(n,sz) R_alloc( n, sz )
#define nfree(...) nofun()
#define check_interrupt R_CheckUserInterrupt

#else

#include <stdio.h>
#include <stdlib.h>

#define mess(...) fprintf(stderr,__VA_ARGS__)

#define nalloc(n,sz) malloc( (size_t)(sz) * (n)) 
#define nfree(p) free(p)
#define check_interrupt(...) nofun()

#endif // NOT_R_PACKAGE

#endif // _ETC_H_
