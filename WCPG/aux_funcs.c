#include "aux_funcs.h"

void *(*wcpg_malloc)(size_t) = NULL;
void *(*wcpg_calloc)(size_t, size_t) = NULL;
void *(*wcpg_realloc)(void *, size_t) = NULL;
void (*wcpg_free)(void *) = NULL;

static inline void *wcpg_effective_malloc(size_t size) {
  if (wcpg_malloc == NULL) return malloc(size);
  return wcpg_malloc(size);
}

static inline void *wcpg_effective_calloc(size_t nmemb, size_t size) {
  if (wcpg_calloc == NULL) return calloc(nmemb, size);
  return wcpg_calloc(nmemb, size);
}

static inline void *wcpg_effective_realloc(void *ptr, size_t size) {
  if (wcpg_realloc == NULL) return realloc(ptr, size);
  return wcpg_realloc(ptr, size);
}

static inline void wcpg_effective_free(void *ptr) {
  if (wcpg_free == NULL) return free(ptr);
  return wcpg_free(ptr);
}

void wcpg_set_memory_funcs(void *(*my_malloc)(size_t),
			   void *(*my_calloc)(size_t, size_t),
			   void *(*my_realloc)(void *, size_t),
			   void (*my_free)(void *)) {
  wcpg_malloc = my_malloc;
  wcpg_calloc = my_calloc;
  wcpg_realloc = my_realloc;
  wcpg_free = my_free;
}

void wcpg_get_memory_funcs(void *(**a_malloc)(size_t),
			   void *(**a_calloc)(size_t, size_t),
			   void *(**a_realloc)(void *, size_t),
			   void (**a_free)(void *)) {
  *a_malloc = wcpg_malloc;
  *a_calloc = wcpg_calloc;
  *a_realloc = wcpg_realloc;
  *a_free = wcpg_free;
}

/* Allocates size bytes of memory; on failure prints an error message
   and calls exit(1)
*/
void *wcpgSafeMalloc(size_t size) {
  void *ptr;

  ptr = wcpg_effective_malloc(size);
  if (ptr == NULL) {
    fprintf(stderr, "Could not allocate %zd bytes of memory\n", size);
    exit(1);
  }

  return ptr;
}

/* Allocates an array of nmemb elements, each of which occupies size
   bytes of memory; on failure prints an error message and calls
   exit(1)
*/
void *wcpgSafeCalloc(size_t nmemb, size_t size)
{
  void *ptr;

  ptr = wcpg_effective_calloc(nmemb, size);
  if (ptr == NULL) {
    fprintf(stderr, "Could not allocate an array of %zd elements with %zd bytes of memory each\n", nmemb, size);
    exit(1);
  }

  return ptr;
}

/* Reallocates size bytes of memory at the location pointed by ptr;
   on failure prints an error message and calls exit(1)
*/
void *wcpgSafeRealloc(void *ptr, size_t size)
{
  void *newPtr;

  newPtr = wcpg_effective_realloc(ptr, size);
  if ((newPtr == NULL) && (!((size == 0) && (ptr != NULL)))) {
      fprintf(stderr, "Could not rellocate %zd bytes of memory at address %p\n", size, ptr);
    exit(1);
  }

  return newPtr;
}

/* Frees the memory at the location pointed by ptr */
void wcpgSafeFree(void *ptr) {
  wcpg_effective_free(ptr);
}
