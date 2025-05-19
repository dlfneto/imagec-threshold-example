#include "alloc.h"

alloc_stack * stack_top = NULL;

void mm_stack_push() {
    alloc_stack *new_stack = (alloc_stack *)malloc(sizeof(alloc_stack));
    if (!new_stack) {
        fprintf(stderr, "Memory allocation failed for stack context.\n");
        exit(1);
    }
    new_stack->head = NULL;
    new_stack->prev = stack_top;
    new_stack->is_stopped = 0;
    stack_top = new_stack;
}

void mm_stack_pop() {
    if (!stack_top) return;
    
    alloc_node *node = stack_top->head;
    while (node) {
        alloc_node *next = node->next;
        free(node->ptr);
        free(node);
        node = next;
    }

    alloc_stack *old_top = stack_top;
    stack_top = stack_top->prev;
    free(old_top);
}

void mm_stack_stop() {
    if (stack_top) {
        stack_top->is_stopped = 1;
    }
}

void mm_stack_resume() {
    if (stack_top) {
        stack_top->is_stopped = 0;
    }
}

void mm_register_alloc(void *ptr) {
    if (!ptr || !stack_top || stack_top->is_stopped) return;

    alloc_node *new_node = (alloc_node *)malloc(sizeof(alloc_node));
    if (!new_node) {
        fprintf(stderr, "Memory allocation failed for tracking.\n");
        exit(1);
    }

    new_node->ptr = ptr;
    new_node->next = stack_top->head;
    stack_top->head = new_node;
}

void mm_unregister_alloc(void *ptr) {
    if (!ptr || !stack_top) return;

    alloc_node **current = &(stack_top->head);
    while (*current) {
        if ((*current)->ptr == ptr) {
            alloc_node *to_free = *current;
            *current = (*current)->next;
            free(to_free);
            return;
        }
        current = &((*current)->next);
    }
}

void * mm_alloc_memory(size_t size) {
    void *ptr = malloc(size);
    mm_register_alloc(ptr);
    return ptr;
}

void * mm_realloc_memory(void *ptr, size_t size) {
    void *new_ptr = realloc(ptr, size);
    if (new_ptr) {
        mm_unregister_alloc(ptr);
        mm_register_alloc(new_ptr);
    }
    return new_ptr;
}
