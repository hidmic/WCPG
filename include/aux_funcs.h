/*
 * aux_funcs.h
 *
 *This file contains headers for the safe memory allocation functions
 *common to all project source files.
 */

#ifndef WCPG_AUX_FUNCS_H_
#define WCPG_AUX_FUNCS_H_



#ifdef __cplusplus
extern "C" {
#endif


#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <math.h>

void *wcpgSafeMalloc(size_t size);
void *wcpgSafeCalloc(size_t nmemb, size_t size);
void *wcpgSafeRealloc(void *ptr, size_t size);
void wcpgSafeFree(void *ptr);




#ifdef __cplusplus
}
#endif






#endif /* WCPG_AUX_FUNCS_H_ */
