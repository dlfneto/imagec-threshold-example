#ifndef __ARRAY_H__
#define __ARRAY_H__

/**
  struct dynamic_array {
    int count;
    int capacity;
    void * items;
  };
*/

typedef struct double_array {
  int count;
  int capacity;
  double * items;
} double_arr;
#define DOUBLE_ARR(size) (double_arr) {.count = 0, .capacity = (size), .items=ALLOC(sizeof(double)*(size)) }

typedef struct integer_array {
  int count;
  int capacity;
  double * items;
} int_arr;
#define INT_ARR(size) (int_arr) {.count = 0, .capacity = (size), .items=ALLOC(sizeof(int)*(size)) }


#define VGET(arr, idx) arr.items[idx]

#define APPEND(arr, item) do {\
    if((arr).capacity <= (arr).count + 1){\
        int new_cap = (arr).capacity * 1.5 + 7;\
        (arr).items = REALLOC((arr).items, sizeof(*(arr).items)*new_cap);\
        assert((arr).items);\
        (arr).capacity = new_cap; \
    }\
    (arr).items[(arr).count++] = item;\
} while(0)

#define DEQUEUE(arr, dest) do {\
  dest = (arr).items[(arr).count-1];\
  (arr).count--;\
} while(0)

#define PUSH(arr, item) do {\
  if((arr).capacity <= (arr).count + 1){\
      int new_cap = (arr).capacity * 1.5 + 7;\
      (arr).items = REALLOC((arr).items, sizeof(*(arr).items)*new_cap);\
      assert((arr).items);\
      (arr).capacity = new_cap; \
  }\
  for(int i = (arr).count - 1; i >= 0; --i){\
    (arr).items[i + 1] = (arr).items[i];\
  }\
  (arr).items[0] = (item);\
  (arr).count++;\
} while(0)

#define POP(arr, dest) do {\
  assert((arr).count > 0);\
  dest = (arr).items[0]; \
  for(int i = 0; i < (arr).count - 1; ++i){\
    (arr).items[i] = (arr).items[i + 1];\
  }\
  (arr).count--; \
} while(0)

#endif /* __ARRAY_H__ */
