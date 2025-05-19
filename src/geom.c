#include "geom.h"

void print_polygon(polygon_t p){
    for(int i = 0; i < p.count; ++i){
        printf("(%d, %d)\n", p.points.items[i].x, p.points.items[i].y);
    }
}

vec2_t vec2(int x, int y){
    return (vec2_t) {.x=x, .y=y};
}

geometry_t rect(int x, int y, int x2, int y2){
    rect_t r = {
        .x1 = x,
        .x2 = x2,
        .y1 = y,
        .y2 = y2
    };
    return (geometry_t) {.type = G_RECT, .rect = r,};
}

geometry_t polygon(int count, ...){
    assert(count > 2);
    va_list args;
    va_start(args, count);
    polygon_t poly = {0};
    poly.count = count;
    poly.points = VEC2_ARR(count);
    for(int i = 0; i < count; ++i){
        vec2_t item = va_arg(args, vec2_t);
        APPEND(poly.points, item);
    }
    va_end(args);
    return (geometry_t) {.type = G_POLYGON, .polygon = poly};
}


int polygon_is_point_inside(polygon_t polygon, vec2_t p) {
    int i, j;
    int inside = 0;
    int n = polygon.count;
    vec2_t * points = polygon.points.items;
    int intersect;

    int y = p.y;
    int x = p.x;

    for (i = 0, j = n - 1; i < n; j = i++) {
        int xi = points[i].x, yi = points[i].y;
        int xj = points[j].x, yj = points[j].y;
        intersect = ((yi > y) != (yj > y)) &&
                         (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
        if (intersect) {
            inside = !inside;
        }
    }

    return inside;
}