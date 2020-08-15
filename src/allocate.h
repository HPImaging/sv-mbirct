
#ifndef ALLOCATE_INC
#define ALLOCATE_INC

#include <stdlib.h>

void *get_spc(size_t num, size_t size);
void *mget_spc(size_t num, size_t size);
void **get_img(int wd,int ht, size_t size);
void ***get_3D(int N, int M, int A, size_t size);
void free_img(void **pt);
void free_3D(void ***pt);
void *multialloc(size_t s, int d, ...);
void multifree(void *r,int d);


#endif
