#ifndef __GEOM_IMAGEC_H__
#define __GEOM_IMAGEC_H__

#include <stdarg.h>
#include <assert.h>
#include "alloc.h"
#include "array.h"

#define VEC2_ARR(size) (vec2_arr) {.count = 0, .capacity = (size), .items=ALLOC(sizeof(vec2_arr)*(size)) }

typedef enum geom_type {G_RECT, G_POLYGON} geom_type;
typedef struct bounding_rect {
    int x1;
    int x2;
    int y1;
    int y2;
} rect_t;

typedef struct vec2 {
    int x;
    int y;
} vec2_t;

typedef struct vec2_array {
    int count;
    int capacity;
    vec2_t * items;
} vec2_arr;

typedef struct bounding_polygon {
    int count;
    vec2_arr points;
} polygon_t;

typedef struct geometry {
    geom_type type;
    union {
        rect_t rect;
        polygon_t polygon;
    };
} geometry_t;


void print_polygon(polygon_t p);
vec2_t vec2(int x, int y);
geometry_t polygon(int count, ...);
geometry_t rect(int x, int y, int x2, int y2);
int polygon_is_point_inside(polygon_t polygon, vec2_t p);

#endif /* __GEOM_IMAGEC_H__ */