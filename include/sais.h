/*
 * sais.h for sais
 * Copyright (c) 2008-2010 Yuta Mori All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _SAIS_H
#define _SAIS_H 1

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <inttypes.h>

#ifndef SAIS_API
# define SAIS_API 
#endif

/*- Datatypes -*/
#ifndef _SA_UINT8_T
#define _SA_UINT8_T
typedef uint8_t sa_uint8_t;
#endif
#ifndef _SA_UINT16_T
#define _SA_UINT16_T
typedef uint16_t sa_uint16_t;
#endif
#ifndef _SA_INT32_T
#define _SA_INT32_T
typedef int32_t sa_int32_t;
#endif
#ifndef SA_PRIdINT32
#define SA_PRIdINT32 PRId32
#endif


/*- Prototypes -*/

/**
 * Constructs the suffix array of a given string.
 * @param T[0..n-1] The input string.
 * @param SA[0..n-1] The output array of suffixes.
 * @param n The length of the given string.
 * @param k The alphabet size.
 * @return 0 if no error occurred, -1 or -2 otherwise.
 */
SAIS_API
sa_int32_t
sais_u8(const sa_uint8_t *T, sa_int32_t *SA, sa_int32_t n, sa_int32_t k);

SAIS_API
sa_int32_t
sais_u16(const sa_uint16_t *T, sa_int32_t *SA, sa_int32_t n, sa_int32_t k);

SAIS_API
sa_int32_t
sais_i32(const sa_int32_t *T, sa_int32_t *SA, sa_int32_t n, sa_int32_t k);


/**
 * Constructs the burrows-wheeler transformed string of a given string.
 * @param T[0..n-1] The input string.
 * @param U[0..n-1] The output string. (can be T)
 * @param A[0..n-1] The temporary array. (can be NULL)
 * @param n The length of the given string.
 * @return The primary index if no error occurred, -1 or -2 otherwise.
 */
SAIS_API
sa_int32_t
sais_u8_bwt(const sa_uint8_t *T, sa_uint8_t *U, sa_int32_t *A, sa_int32_t n, sa_int32_t k);

SAIS_API
sa_int32_t
sais_u16_bwt(const sa_uint16_t *T, sa_uint16_t *U, sa_int32_t *A, sa_int32_t n, sa_int32_t k);

SAIS_API
sa_int32_t
sais_i32_bwt(const sa_int32_t *T, sa_int32_t *U, sa_int32_t *A, sa_int32_t n, sa_int32_t k);


/**
 * Returns the version of the sais library.
 * @return The version number string.
 */
SAIS_API
const char *
sais_version(void);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* _SAIS_H */
