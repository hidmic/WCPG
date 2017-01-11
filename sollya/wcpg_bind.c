/*

  Copyright 2007-2016 by

  Laboratoire de l'Informatique du Parallelisme,
  UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668,

  and by

  Laboratoire d'Informatique de Paris 6, equipe PEQUAN,
  UPMC Universite Paris 06 - CNRS - UMR 7606 - LIP6, Paris, France.

  Contributor Ch. Lauter

  christoph.lauter@ens-lyon.org

  This software is a computer program whose purpose is to provide an
  environment for safe floating-point code development. It is
  particularly targeted to the automated implementation of
  mathematical floating-point libraries (libm). Amongst other features,
  it offers a certified infinity norm, an automatic polynomial
  implementer and a fast Remez algorithm.

  This software is governed by the CeCILL-C license under French law and
  abiding by the rules of distribution of free software.  You can  use,
  modify and/ or redistribute the software under the terms of the CeCILL-C
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and,  more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL-C license and that you accept its terms.

  This program is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

*/

#include <mpfr.h>
#include <sollya.h>
#include <wcpg.h>

/* Example of an external procedure linked to an identifier in sollya

   Compile with

   gcc -g -O0 -I/home/lauter/sollya-install/include -I/home/lauter/2016-WCPG/install/include -fPIC -Wall -c wcpg_bind.c && gcc -g -O0 -L/home/lauter/2016-WCPG/install/lib/ -L/home/lauter/sollya-install/lib/ -fPIC -shared -o wcpg_bind.so wcpg_bind.o -lwcpg -lsollya

   Procedure __wcpg will be linked by

   externalproc(__wcpg, "wcpg_bind.so", (list of constant, list of constant, list of constant, list of constant, integer, integer, integer, constant) -> list of constant);

*/

int __wcpg_wrapper(mpfr_t *W, mpfr_t *A, mpfr_t *B, mpfr_t *C, mpfr_t *D, int n, int p, int q, mpfr_t eps) {
  return !!WCPG_mp(W, A, B, C, D, (uint64_t) n, (uint64_t) p, (uint64_t) q, eps);
}

mp_prec_t getSollyaPrec() {
  sollya_obj_t t;
  int64_t p;
  
  t = sollya_lib_get_prec();
  if (!sollya_lib_get_constant_as_int64(&p, t)) {
    p = 12;
  }
  sollya_lib_clear_obj(t);
  return (mp_prec_t) p;
}

int __wcpg_sollya(sollya_constant_list_t *Wptr,
		  sollya_constant_list_t A,
		  sollya_constant_list_t B,
		  sollya_constant_list_t C,
		  sollya_constant_list_t D,
		  int n, int p, int q, mpfr_t eps) {
  int okay;
  mpfr_t *A_mat, *B_mat, *C_mat, *D_mat, *W_mat, *temp;
  sollya_constant_list_t W, curr;
  mp_prec_t prec;
  int i, j;

  prec = getSollyaPrec();

  A_mat = (mpfr_t *) sollya_lib_calloc(n * n, sizeof(mpfr_t));
  B_mat = (mpfr_t *) sollya_lib_calloc(n * q, sizeof(mpfr_t));
  C_mat = (mpfr_t *) sollya_lib_calloc(p * n, sizeof(mpfr_t));
  D_mat = (mpfr_t *) sollya_lib_calloc(p * q, sizeof(mpfr_t));
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      mpfr_init2(A_mat[i * n + j], prec);
    }
  }
  for (i=0;i<n;i++) {
    for (j=0;j<q;j++) {
      mpfr_init2(B_mat[i * q + j], prec);
    }
  }
  for (i=0;i<p;i++) {
    for (j=0;j<n;j++) {
      mpfr_init2(C_mat[i * n + j], prec);
    }
  }
  for (i=0;i<p;i++) {
    for (j=0;j<q;j++) {
      mpfr_init2(D_mat[i * q + j], prec);
    }
  }

  for (curr=A,i=0;!sollya_lib_is_empty_constant_list(curr);curr=sollya_lib_get_constant_list_tail(curr),i++) {
    mpfr_set_prec(A_mat[i], mpfr_get_prec(*(sollya_lib_get_constant_list_head(curr))));
    mpfr_set(A_mat[i], *(sollya_lib_get_constant_list_head(curr)), GMP_RNDN); /* exact */
  }
  for (curr=B,i=0;!sollya_lib_is_empty_constant_list(curr);curr=sollya_lib_get_constant_list_tail(curr),i++) {
    mpfr_set_prec(B_mat[i], mpfr_get_prec(*(sollya_lib_get_constant_list_head(curr))));
    mpfr_set(B_mat[i], *(sollya_lib_get_constant_list_head(curr)), GMP_RNDN); /* exact */
  }
  for (curr=C,i=0;!sollya_lib_is_empty_constant_list(curr);curr=sollya_lib_get_constant_list_tail(curr),i++) {
    mpfr_set_prec(C_mat[i], mpfr_get_prec(*(sollya_lib_get_constant_list_head(curr))));
    mpfr_set(C_mat[i], *(sollya_lib_get_constant_list_head(curr)), GMP_RNDN); /* exact */
  }
  for (curr=D,i=0;!sollya_lib_is_empty_constant_list(curr);curr=sollya_lib_get_constant_list_tail(curr),i++) {
    mpfr_set_prec(D_mat[i], mpfr_get_prec(*(sollya_lib_get_constant_list_head(curr))));
    mpfr_set(D_mat[i], *(sollya_lib_get_constant_list_head(curr)), GMP_RNDN); /* exact */
  }

  W_mat = (mpfr_t *) sollya_lib_calloc(p * q, sizeof(mpfr_t));
  for (i=0;i<p;i++) {
    for (j=0;j<q;j++) {
      mpfr_init2(W_mat[i * q + j], prec);
    }
  }

  okay = __wcpg_wrapper(W_mat, A_mat, B_mat, C_mat, D_mat, n, p, q, eps);

  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      mpfr_clear(A_mat[i * n + j]);
    }
  }
  for (i=0;i<n;i++) {
    for (j=0;j<q;j++) {
      mpfr_clear(B_mat[i * q + j]);
    }
  }
  for (i=0;i<p;i++) {
    for (j=0;j<n;j++) {
      mpfr_clear(C_mat[i * n + j]);
    }
  }
  for (i=0;i<p;i++) {
    for (j=0;j<q;j++) {
      mpfr_clear(D_mat[i * q + j]);
    }
  }  
  sollya_lib_free(D_mat);
  sollya_lib_free(C_mat);
  sollya_lib_free(B_mat);
  sollya_lib_free(A_mat);

  W = NULL;
  for (i=0;i<p;i++) {
    for (j=0;j<q;j++) {
      temp = (mpfr_t *) sollya_lib_malloc(sizeof(mpfr_t));
      mpfr_init2(*temp, mpfr_get_prec(W_mat[i * q + j]));
      mpfr_set(*temp, W_mat[i * q + j], GMP_RNDN); /* exact */
      W = sollya_lib_construct_constant_list(temp, W);
    }
  }  
    
  for (i=0;i<p;i++) {
    for (j=0;j<q;j++) {
      mpfr_clear(W_mat[i * q + j]);
    }
  }  
  sollya_lib_free(W_mat);
  if (okay) {
    *Wptr = W;
  } else {
    sollya_lib_clear_constant_list(W);
  }

  return okay;
}

/* Signature (list of constant, list of constant, list of constant, list of constant, integer, integer, integer, constant) -> list of constant */
int __wcpg(sollya_constant_list_t *res, void **args) {
  int okay;
  void *(*old_malloc)(size_t);
  void *(*old_calloc)(size_t, size_t);
  void *(*old_realloc)(void *, size_t);
  void (*old_free)(void *);
  
  okay = 0;
  wcpg_get_memory_funcs(&old_malloc, &old_calloc, &old_realloc, &old_free);
  wcpg_set_memory_funcs(sollya_lib_malloc, sollya_lib_calloc, sollya_lib_realloc, sollya_lib_free);

  okay = __wcpg_sollya((sollya_constant_list_t *) res,
		       ((sollya_constant_list_t) (args[0])),
		       ((sollya_constant_list_t) (args[1])),
		       ((sollya_constant_list_t) (args[2])),
		       ((sollya_constant_list_t) (args[3])),
		       *((int *) (args[4])),
		       *((int *) (args[5])),
		       *((int *) (args[6])),
		       *((mpfr_t *) (args[7])));
    
  wcpg_set_memory_funcs(old_malloc, old_calloc, old_realloc, old_free);
  return okay;
}
