#ifndef __ALLOC_H__
#define __ALLOC_H__

#include <stdlib.h>
#include <stdio.h>

typedef struct alloc_node {
    void *ptr;
    struct alloc_node *next;
} alloc_node;

typedef struct alloc_stack {
    alloc_node *head;
    struct alloc_stack *prev;
    int is_stopped;
} alloc_stack;

void mm_stack_push();
void mm_stack_pop();
void mm_stack_stop();
void mm_stack_resume();
void * mm_alloc_memory(size_t size);
void * mm_realloc_memory(void *ptr, size_t size);
void mm_register_alloc(void *ptr);
void mm_unregister_alloc(void *ptr);

/**
 * Create new MM_BEGIN scope & Begin tracking next memory allocations.
 */
#define MM_BEGIN() mm_stack_push()

/**
 * Clean current tracked memory in MM_BEGIN scope.
 */
#define MM_CLEAN() mm_stack_pop()

/**
 * Stop tracking allocation of current MM_BEGIN scope.
 */
#define MM_STOP() mm_stack_stop()

/**
 * Resume tracking allocation of current MM_BEGIN scope.
 */
#define MM_RESUME() mm_stack_resume()

#define ALLOC(size) mm_alloc_memory(size)
#define REALLOC(ptr, size) mm_realloc_memory(ptr, size)

#define FREE(ptr) ({ \
    if (ptr) { \
        mm_unregister_alloc(ptr); \
        free(ptr); \
        ptr = NULL; \
    } \
})

#endif /* __ALLOC_H__ */
