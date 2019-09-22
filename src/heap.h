#include <stdio.h>
#include <stdlib.h>

struct heap {
  int size;
  int alloc;
  struct heap_node{
    float x;
    int pt;  /* may be used as array pointer or data */
  }*node;
};



void initialize_heap(struct heap *h);
void free_heap(struct heap *h);
void heap_insert(struct heap *h, struct heap_node *p);
void get_heap_max(struct heap *h, struct heap_node *p);
void maintain_heap(struct heap *h, int n);
