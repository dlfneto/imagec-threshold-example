
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "imagec.h"

typedef unsigned char u8_t;
typedef struct pixel {
    union
    {
        u8_t lum;
        struct {
            u8_t r;
            u8_t g;
            u8_t b;
            u8_t a;
        };
    };
} pixel_t;

pixel_t ** data_to_pixel(int w, int h, void * data, int channels){
    pixel_t ** pixels = (pixel_t **) ALLOC(sizeof(*pixels)*h);
    assert(pixels);
    for(int i = 0; i < h; ++i){
        pixels[i] = (pixel_t *) ALLOC(sizeof(**pixels)*w);
        assert(pixels[i]);
    }
    size_t offset_pxl = sizeof(u8_t);
    int c = 0;
    int c_offset = 0;
    void * pixel;
    for(int i = 0; i < h; ++i){
        for(int j = 0; j < w; ++j){
            c = channels;
            while(c > 0){
                pixel = &pixels[i][j];
                c_offset = (channels-c);
                *(u8_t *)(pixel+(offset_pxl*c_offset)) = *(u8_t *)(data+((i*w + j)*channels+c_offset));
                c -= 1;
            }
        }
    }
    return pixels;
}

image_t create_image(int w, int h, int ch){
    image_t image = {
        .h = h,
        .w = w,
        .channels = ch
    };
    image.data = ALLOC(sizeof(*image.data)*image.w*image.h*image.channels);
    assert(image.data);
    return image;
}

image_t load_image(char * filename){
    image_t image = {0};
    void * data = stbi_load(filename, &image.w, &image.h, &image.channels, 0);
    image.data = ALLOC(sizeof(*image.data)*image.w*image.h*image.channels);
    assert(image.data);

    int size = image.w*image.h*image.channels;
    unsigned char * value;
    for(int i = 0; i < size; ++i){
        value = data+i;
        image.data[i] = *value;
    }
    FREE(data);
    return image;
}

void write_image(char * filename, image_t image){
    int size = image.w*image.h*image.channels;
    unsigned char * data = ALLOC(sizeof(*data)*size);
    assert(data);

    long double value = 0;
    for(int i = 0; i < size; ++i){
        value = image.data[i];
        if(value > 255) value = 255;
        if(value < 0) value = 0;
        data[i] = value;
    }
    stbi_write_png(filename, image.w, image.h, image.channels, data, 0);
    FREE(data);
}

void free_image(image_t * image){
    image->h = 0;
    image->w = 0;
    image->channels = 0;
    if(image->data) FREE(image->data);
    image->data = NULL;
}

image_t clone(image_t source){
    image_t new_image = {
        .channels = source.channels,
        .h = source.h,
        .w = source.w,
        .data = ALLOC(sizeof(*source.data)*source.w*source.h*source.channels)
    };
    assert(new_image.data);
    int size = source.w*source.h*source.channels;
    for(int i = 0; i < size; ++i) new_image.data[i] = source.data[i];
    return new_image;
}

image_t rgb2gray(image_t image){
    if(image.channels < 3) return (image_t){0};
    image_t new_image = {
        .h = image.h,
        .w = image.w,
        .data = NULL,
        .channels = 1
    };

    new_image.data = ALLOC(sizeof(*new_image.data)*new_image.h*new_image.w);
    assert(new_image.data);

    size_t max = new_image.h*new_image.w;
    for(size_t i = 0; i < max; ++i){
        if(image.data[(i*image.channels)] == VOID_PIXEL){
            new_image.data[i] = VOID_PIXEL;
            continue;
        }
        new_image.data[i] = 0.3 * image.data[(i*image.channels)] + 0.59 * image.data[(i*image.channels)+1] + 0.11 * image.data[(i*image.channels)+2];
        // new_image.data[i] = (image.data[(i*image.channels)] + image.data[(i*image.channels)+1] + image.data[(i*image.channels)+2])/3;
    }

    return new_image;
}

image_t crop_rect(image_t image, rect_t bounds);
image_t crop_polygon(image_t image, polygon_t bounds);
image_t cut_rect(image_t image, rect_t bounds);
image_t cut_polygon(image_t image, polygon_t bounds);

image_t crop(image_t image, geometry_t geometry_bounds){
    if(geometry_bounds.type == G_RECT){
        rect_t r = geometry_bounds.rect;
        return crop_rect(image, r);
    }
    if(geometry_bounds.type == G_POLYGON){
        polygon_t p = geometry_bounds.polygon;
        return crop_polygon(image, p);
    }
    return (image_t) {0};
}

image_t cut(image_t image, geometry_t geometry_bounds){
    if(geometry_bounds.type == G_RECT){
        rect_t r = geometry_bounds.rect;
        return cut_rect(image, r);
    }
    if(geometry_bounds.type == G_POLYGON){
        polygon_t p = geometry_bounds.polygon;
        return cut_polygon(image, p);
    }
    return (image_t) {0};
}

image_t cut_rect(image_t image, rect_t bounds){
    int swap = 0;

    int x1, x2, y1, y2;
    x1 = bounds.x1;
    x2 = bounds.x2;
    y1 = bounds.y1;
    y2 = bounds.y2;

    if(x1 > x2){
        swap = x1;
        x1 = x2;
        x2 = swap;
    }

    if(y1 > y2){
        swap = y1;
        y1 = y2;
        y2 = swap;
    }

    int h = y2-y1;
    int w = x2-x1;

    image_t cropped = {
        .h = h,
        .w = w,
        .channels = image.channels,
        .data = NULL
    };
    cropped.data = ALLOC(sizeof(*cropped.data)*h*w*cropped.channels);

    size_t pos = 0;
    int tx = 0, ty = 0;

    for(size_t i = 0; i < h; ++i){
        ty = (y1 + i);
        for(size_t j = 0; j < w; ++j){
            tx = (x1 + j);
            for(size_t c = 0; c < image.channels; ++c){
                pos = (ty * image.w + tx) * image.channels + c;
                if(pos > (image.w * image.h * image.channels)){
                    PXL_AT(cropped, j, i, c) = 0;
                }
                else {
                    PXL_AT(cropped, j, i, c) = image.data[pos];
                    image.data[pos] = VOID_PIXEL;
                }
            }
        }
    }
    return cropped;
}

image_t crop_rect(image_t image, rect_t bounds){
    int swap = 0;

    int x1, x2, y1, y2;
    x1 = bounds.x1;
    x2 = bounds.x2;
    y1 = bounds.y1;
    y2 = bounds.y2;

    if(x1 > x2){
        swap = x1;
        x1 = x2;
        x2 = swap;
    }

    if(y1 > y2){
        swap = y1;
        y1 = y2;
        y2 = swap;
    }

    if (x1 < 0) x1 = 0;
    if (x2 >= image.w) x2 = image.w - 1;
    if (y1 < 0) y1 = 0;
    if (y2 >= image.h) y2 = image.h - 1;

    int h = y2 - y1 + 1;
    int w = x2 - x1 + 1;

    image_t cropped = {
        .h = h,
        .w = w,
        .channels = image.channels,
        .data = NULL
    };
    cropped.data = ALLOC(sizeof(*cropped.data)*h*w*cropped.channels);

    size_t pos = 0;
    int tx = 0, ty = 0;

    for(size_t i = 0; i < h; ++i){
        ty = (y1 + i);
        for(size_t j = 0; j < w; ++j){
            tx = (x1 + j);
            for(size_t c = 0; c < image.channels; ++c){
                pos = (ty * image.w + tx) * image.channels + c;
                if(pos >= (image.w * image.h * image.channels))
                    PXL_AT(cropped, j, i, c) = 0;
                else
                    PXL_AT(cropped, j, i, c) = image.data[pos];
            }
        }
    }
    return cropped;
}

image_t cut_polygon(image_t image, polygon_t bounds){
    int max_w = 0;
    int min_w = INT_MAX;
    int max_h = 0;
    int min_h = INT_MAX;

    for(int i = 0; i < bounds.count; ++i){
        if(bounds.points.items[i].x > max_w) max_w = bounds.points.items[i].x;
        if(bounds.points.items[i].x < min_w) min_w = bounds.points.items[i].x;
        if(bounds.points.items[i].y > max_h) max_h = bounds.points.items[i].y;
        if(bounds.points.items[i].y < min_h) min_h = bounds.points.items[i].y;
    }

    int h = max_h - min_h;
    int w = max_w - min_w;

    image_t cropped = {
        .h = h,
        .w = w,
        .channels = image.channels,
        .data = NULL
    };

    cropped.data = ALLOC(sizeof(*cropped.data)*max_w*max_h*cropped.channels);
    size_t pos = 0;
    int tx = 0, ty = 0;

    for(size_t i = 0; i < h; ++i){
        ty = (min_h + i);
        for(size_t j = 0; j < w; ++j){
            tx = (min_w + j);
            if(!polygon_is_point_inside(bounds, vec2(tx, ty))){
                for(size_t c = 0; c < image.channels; ++c) PXL_AT(cropped, j, i, c) = VOID_PIXEL;
                if(image.channels == 4) PXL_AT(cropped, j, i, 3) = 0;
                continue;
            }
            for(size_t c = 0; c < image.channels; ++c){
                pos = (ty * image.w + tx) * image.channels + c;
                if(pos > (image.w * image.h * image.channels)){
                    PXL_AT(cropped, j, i, c) = 0;
                }
                else {
                    PXL_AT(cropped, j, i, c) = image.data[pos];
                    image.data[pos] = VOID_PIXEL;
                }
            }
        }
    }
    return cropped;
}

image_t crop_polygon(image_t image, polygon_t bounds){
    int max_w = 0;
    int min_w = INT_MAX;
    int max_h = 0;
    int min_h = INT_MAX;

    for(int i = 0; i < bounds.count; ++i){
        if(bounds.points.items[i].x > max_w) max_w = bounds.points.items[i].x;
        if(bounds.points.items[i].x < min_w) min_w = bounds.points.items[i].x;
        if(bounds.points.items[i].y > max_h) max_h = bounds.points.items[i].y;
        if(bounds.points.items[i].y < min_h) min_h = bounds.points.items[i].y;
    }

    int h = max_h - min_h;
    int w = max_w - min_w;

    image_t cropped = {
        .h = h,
        .w = w,
        .channels = image.channels,
        .data = NULL
    };

    cropped.data = ALLOC(sizeof(*cropped.data)*max_w*max_h*cropped.channels);
    size_t pos = 0;
    int tx = 0, ty = 0;

    for(size_t i = 0; i < h; ++i){
        ty = (min_h + i);
        for(size_t j = 0; j < w; ++j){
            tx = (min_w + j);
            if(!polygon_is_point_inside(bounds, vec2(tx, ty))){
                for(size_t c = 0; c < image.channels; ++c) PXL_AT(cropped, j, i, c) = VOID_PIXEL;
                if(image.channels == 4) PXL_AT(cropped, j, i, 3) = 0;
                continue;
            }
            for(size_t c = 0; c < image.channels; ++c){
                pos = (ty * image.w + tx) * image.channels + c;
                if(pos > (image.w * image.h * image.channels))
                    PXL_AT(cropped, j, i, c) = 0;
                else
                    PXL_AT(cropped, j, i, c) = image.data[pos];
            }
        }
    }
    return cropped;
}

void upgrade(image_t * image, int channels){
    if(channels <= image->channels) return;
    image_t clonned = clone(*image);

    image->w = clonned.w;
    image->h = clonned.h;
    image->channels = channels;
    image->data = ALLOC(sizeof(*image->data)*clonned.h*clonned.w*channels);
    assert(image->data);

    FOREACH_PXL(*image, {
        if(c < clonned.channels)
            PXL_AT(*image, x, y, c) = PXL_AT(clonned, x, y, c);
        else
            PXL_AT(*image, x, y, c) = 0;
    });

    free_image(&clonned);
}

image_t paste(image_t src, image_t dest, int x, int y){

    if(dest.channels < src.channels) upgrade(&dest, src.channels);

    size_t pos = 0;
    int tx = 0, ty = 0;

    int channels = dest.channels;
    
    for(size_t i = 0; i < src.h; ++i){
        ty = (y + i);
        for(size_t j = 0; j < src.w; ++j){
            tx = (x + j);
            for(size_t c = 0; c < channels; ++c){
                if(src.channels <= c) {
                    if(src.data[(i*src.w + j)*src.channels + 0] == VOID_PIXEL) continue; 
                    PXL_AT(dest, tx, ty, c) = src.data[(i*src.w + j)*src.channels + 0];
                    continue;
                }
                if(src.data[(i*src.w + j)*src.channels + c] == VOID_PIXEL) {
                    continue;
                }
                pos = (ty * dest.w + tx) * dest.channels + c;
                if(pos > (dest.w * dest.h * dest.channels))
                    continue;
                else
                    PXL_AT(dest, tx, ty, c) = src.data[(i*src.w + j)*src.channels + c];
            }
        }
    }
    return dest;
}

image_t convolve(image_t image, matrix_t kernel){
    image_t new_image = clone(image);

    int kxstart = -(kernel.c/2);
    int kystart = -(kernel.r/2);
    int kxend = (kernel.c/2);
    int kyend = (kernel.r/2);

    double value = 0;
    int channels = image.channels;
    int px = 0;
    int py = 0;

    for(size_t y = 0; y < new_image.h; ++y){
        for(size_t x = 0; x < new_image.w; ++x){

            for(size_t c = 0; c < channels; ++c){

                // Ignore alpha channel
                if(PXL_AT(new_image, x, y, c) == VOID_PIXEL) continue;
                if(c == 3) {
                    PXL_AT(new_image, x, y, c) = PXL_AT(image, x, y, c);
                    continue;
                }

                value = 0;

                for(int kx = kxstart; kx <= kxend; ++kx){
                    px = x + kx;
                    if (px < 0) px = 0;
                    if (px >= image.w) px = image.w-1;
                    for(int ky = kystart; ky <= kyend; ++ky){

                        py = y + ky;
                        if (py < 0) py = 0;
                        if (py >= image.h) py = image.h-1;   
                        
                        /**
                         * For VOID_PIXELS found, ignore; Then to
                         * avoid side effects multiply by the central
                         * pixel.
                         */
                        if(PXL_AT(image, px, py, c) == VOID_PIXEL){
                            value += PXL_AT(image, x, y, c) * MAT_AT(kernel, kx + kxend, ky + kyend);
                            continue;
                        }

                        value += PXL_AT(image, px, py, c) * MAT_AT(kernel, kx + kxend, ky + kyend);
                    }
                }
                PXL_AT(new_image, x, y, c) = value;
            }
        }
    }

    return new_image;
}

image_t laplacian(image_t image) {
    image_t lapl_image = zeros(image.w, image.h, image.channels);

    matrix_t kernel = laplacian_kernel();
    int kxstart = -(kernel.c / 2);
    int kystart = -(kernel.r / 2);
    int kxend = (kernel.c / 2);
    int kyend = (kernel.r / 2);

    double value = 0.0;
    int py, px;
    for (size_t y = 0; y < image.h; ++y) {
        for (size_t x = 0; x < image.w; ++x) {
            for (size_t c = 0; c < image.channels; ++c) {
                if (c == 3) {
                    PXL_AT(lapl_image, x, y, c) = PXL_AT(image, x, y, c);
                    continue;
                }

                value = 0.0;
                for (int ky = kystart; ky <= kyend; ++ky) {
                    py = y + ky;
                    if (py < 0) py = 0;
                    if (py >= image.h) py = image.h - 1;

                    for (int kx = kxstart; kx <= kxend; ++kx) {
                        px = x + kx;
                        if (px < 0) px = 0;
                        if (px >= image.w) px = image.w - 1;

                        value += PXL_AT(image, px, py, c) * MAT_AT(kernel, kx + kxend, ky + kyend);
                    }
                }
                PXL_AT(lapl_image, x, y, c) = value;
            }
        }
    }

    return lapl_image;
}

image_t zero_crossing(image_t laplacian_img, double threshold) {
    image_t output = zeros(laplacian_img.w, laplacian_img.h, laplacian_img.channels);

    for (size_t y = 1; y < laplacian_img.h - 1; ++y) {
        for (size_t x = 1; x < laplacian_img.w - 1; ++x) {
            for (size_t c = 0; c < laplacian_img.channels; ++c) {
                if (c == 3) {
                    PXL_AT(output, x, y, c) = PXL_AT(laplacian_img, x, y, c);
                    continue;
                }

                double center = PXL_AT(laplacian_img, x, y, c);
                int has_zero_crossing = 0;

                int dx[] = {1, -1, 0, 0};
                int dy[] = {0, 0, 1, -1};

                for (int k = 0; k < 4; ++k) {
                    double neighbor = PXL_AT(laplacian_img, x + dx[k], y + dy[k], c);
                    if ((center < 0 && neighbor > 0) || (center > 0 && neighbor < 0)) {
                        double magnitude = fabs(center - neighbor);
                        if (magnitude > threshold) {
                            has_zero_crossing = 1;
                            break;
                        }
                    }
                }

                PXL_AT(output, x, y, c) = has_zero_crossing ? 255.0 : 0.0;
            }
        }
    }
    return output;
}

image_t laplacian_mapped(image_t image) {
    image_t new_image = laplacian(image);
    for (size_t y = 0; y < new_image.h; ++y) {
        for (size_t x = 0; x < new_image.w; ++x) {
            for (size_t c = 0; c < new_image.channels; ++c) {
                if (c == 3) continue;
                double value = PXL_AT(new_image, x, y, c) + 128.0;
                PXL_AT(new_image, x, y, c) = value;
            }
        }
    }
    return new_image;
}

image_t sobel(image_t image) {
    image_t new_image = clone(image);
    image_t sobel_x = convolve(image, sobel_x_kernel());
    image_t sobel_y = convolve(image, sobel_y_kernel());
    double gx, gy, magnitude;
    for (size_t y = 0; y < new_image.h; ++y) {
        for (size_t x = 0; x < new_image.w; ++x) {
            for (int c = 0; c < image.channels; ++c) {
                // Ignore alpha channel
                if(PXL_AT(new_image, x, y, c) == VOID_PIXEL) continue;
                if (c == 3) {  
                    PXL_AT(new_image, x, y, c) = PXL_AT(image, x, y, c);
                    continue;
                }
                gx = PXL_AT(sobel_x, x, y, c);
                gy = PXL_AT(sobel_y, x, y, c);
                magnitude = sqrt(gx * gx + gy * gy);
                PXL_AT(new_image, x, y, c) = magnitude;
            }
        }
    }
    return new_image;
}

image_t sobel_x5(image_t image) {
    image_t new_image = clone(image);
    image_t sobel_x = convolve(image, sobel_x_kernel_x5());
    image_t sobel_y = convolve(image, sobel_y_kernel_x5());
    for (size_t y = 0; y < new_image.h; ++y) {
        for (size_t x = 0; x < new_image.w; ++x) {
            for (int c = 0; c < image.channels; ++c) {

                if(PXL_AT(new_image, x, y, c) == VOID_PIXEL) continue;

                // Ignore alpha channel
                if (c == 3) {  
                    PXL_AT(new_image, x, y, c) = PXL_AT(image, x, y, c);
                    continue;
                }

                double gx = PXL_AT(sobel_x, x, y, c);
                double gy = PXL_AT(sobel_y, x, y, c);

                double magnitude = sqrt(gx * gx + gy * gy);

                PXL_AT(new_image, x, y, c) = magnitude;
            }
        }
    }
    return new_image;
}

int compare_doubles(const void *a, const void *b) {
    double diff = (*(double *)a) - (*(double *)b);
    return (diff > 0) - (diff < 0);
}


image_t threshold_by(image_t image, double value) {
    /**
     * Threshold based on start index percentile
     *  */ 
    image_t result = clone(image);
    for (size_t y = 0; y < image.h; ++y) {
        for (size_t x = 0; x < image.w; ++x) {
            for (int c = 0; c < image.channels; ++c) {
                if (c == 3) {
                    PXL_AT(result, x, y, c) = PXL_AT(image, x, y, c);
                    continue;
                }
                double val = PXL_AT(image, x, y, c);
                PXL_AT(result, x, y, c) = (val >= value) ? val : 0.0;
            }
        }
    }
    return result;
}

image_t threshold_percentile(image_t gradient_image, double percentile) {
    size_t total_pixels = gradient_image.w * gradient_image.h * gradient_image.channels;
    double * magnitudes = ALLOC(sizeof(double) * total_pixels);
    size_t count = 0;

    /**
     * Create an array with all image data
     *  */ 
    for (size_t y = 0; y < gradient_image.h; ++y) {
        for (size_t x = 0; x < gradient_image.w; ++x) {
            for (int c = 0; c < gradient_image.channels; ++c) {
                // Ignore Alpha Channel
                if (c == 3) continue;

                double val = PXL_AT(gradient_image, x, y, c);

                // Ignore VOID_PIXELS
                if(val == VOID_PIXEL) continue;

                magnitudes[count++] = val;
            }
        }
    }

    /**
     * Sort the array to be able to get the percentile
     *  */ 
    qsort(magnitudes, count, sizeof(double), compare_doubles);

    size_t index = (size_t)(percentile * count);
    if (index >= count) index = count - 1;
    double threshold = magnitudes[index];
    FREE(magnitudes);
    image_t result = clone(gradient_image);

    /**
     * Threshold based on start index percentile
     *  */ 
    for (size_t y = 0; y < result.h; ++y) {
        for (size_t x = 0; x < result.w; ++x) {
            for (int c = 0; c < result.channels; ++c) {
                if (c == 3) {
                    PXL_AT(result, x, y, c) = PXL_AT(gradient_image, x, y, c);
                    continue;
                }
                double val = PXL_AT(gradient_image, x, y, c);
                PXL_AT(result, x, y, c) = (val >= threshold) ? 255.0 : 0.0;
            }
        }
    }

    return result;
}

image_t normalize_nonzero_pixels(image_t input) {
    image_t output = clone(input);
    double min = 1e9, max = -1e9;
    for (int y = 0; y < input.h; ++y) {
        for (int x = 0; x < input.w; ++x) {
            double val = PXL_AT(input, x, y, 0);
            if (val > 0.0) {
                if (val < min) min = val;
                if (val > max) max = val;
            }
        }
    }
    double range = max - min;
    if (range == 0.0) range = 1.0;
    for (int y = 0; y < input.h; ++y) {
        for (int x = 0; x < input.w; ++x) {
            double val = PXL_AT(input, x, y, 0);
            if (val > 0.0) {
                PXL_AT(output, x, y, 0) = 255.0 * (val - min) / range;
            } else {
                PXL_AT(output, x, y, 0) = 0.0;
            }
        }
    }
    return output;
}

int * histdata(image_t image){
    assert(image.channels == 1);
    int * histogram = ALLOC(256 * sizeof(*histogram));
    for(int i = 0; i < 256; ++i) histogram[i] = 0;

    int val = 0;
    for (int y = 0; y < image.h; ++y) {
        for (int x = 0; x < image.w; ++x) {
            for (int c = 0; c < image.channels; ++c) {
                if (c == 3) continue;
                val = (int)PXL_AT(image, x, y, c);
                if (val < 0) val = 0;
                if (val > 255) val = 255;
                histogram[val]++;
            }
        }
    }
    return histogram;
}

int * histdata__ignore_zeros(image_t image){
    assert(image.channels == 1);
    int * histogram = ALLOC(256 * sizeof(*histogram));
    for(int i = 0; i < 256; ++i) histogram[i] = 0.0;
    int val = 0;
    for (int y = 0; y < image.h; ++y) {
        for (int x = 0; x < image.w; ++x) {
            for (int c = 0; c < image.channels; ++c) {
                if (c == 3) continue;
                val = (int) (PXL_AT(image, x, y, c));
                if (val == 0) continue;
                if (val < 0) val = 0;
                if (val > 255) val = 255;
                histogram[val]++;
            }
        }
    }
    return histogram;
}

image_t im2double(image_t image){
    return TRANSFORM(clone(image), { pixel = pixel/255.0 });
}

image_t otsu_from_hist(image_t image, int * histogram){
    
    int total_pixels = 0;
    for(int i = 0; i < 256; ++i) 
        total_pixels += histogram[i];

    double probabilities[256] = {0};
    double sum_all = 0.0;
    for(int i = 0; i < 256; ++i){
        probabilities[i] = (double) (histogram[i]) / (double) total_pixels;
        sum_all += i * probabilities[i];
    }

    double w0 = 0.0;   
    double mu0 = 0.0;  
    double max_sigma = 0.0;
    int best_thresh = 0;
    int consecutive_count = 0;

    for (int t = 0; t < 256; ++t) {
        w0 += probabilities[t];
        mu0 += t * probabilities[t];

        if (w0 == 0.0 || w0 == 1.0)
            continue;

        double w1 = 1.0 - w0;
        double mu1 = (sum_all - mu0) / w1;
        double mu0_norm = mu0 / w0;

        double sigma_b = w0 * w1 * (mu0_norm - mu1) * (mu0_norm - mu1);
        // printf("MAX: %.2f, SIGMA: %.2f | T: %d\n", max_sigma, sigma_b, t);

        if (sigma_b > max_sigma) {
            consecutive_count = 0;
            max_sigma = sigma_b;
            best_thresh = t;
        }

        if (sigma_b == max_sigma) {
            consecutive_count += 1;
            max_sigma = sigma_b;
            best_thresh = t - consecutive_count/2;
        }
    }

    // printf("Best: %d\n", best_thresh);

    image_t result = clone(image);
    for (int y = 0; y < result.h; ++y) {
        for (int x = 0; x < result.w; ++x) {
            for (int c = 0; c < result.channels; ++c) {
                if (c == 3 || PXL_AT(image, x, y, c) == VOID_PIXEL) {
                    PXL_AT(result, x, y, c) = PXL_AT(image, x, y, c);
                    continue;
                }
                int val = (int)PXL_AT(image, x, y, c);
                PXL_AT(result, x, y, c) = (val > best_thresh) ? 255 : 0;
            }
        }
    }

    return result;
}

image_t otsu_threshold(image_t image) {
    int * histogram = histdata(image);
    return otsu_from_hist(image, histogram);
}

void log_image(image_t image){
    printf("W: %d\nH: %d\nCH: %d\n", image.w, image.h, image.channels);
}

image_t resize(image_t src, int w, int h){
    image_t image = create_image(w, h, src.channels);

    double wf = (double) src.w / (double) w;
    double hf = (double) src.h / (double) h;

    int xi, yi;

    for(int y = 0; y < h; ++y){
        for(int x = 0; x < w; ++x){
            xi = (int)(x * wf);
            yi = (int)(y * hf);
            if (xi >= src.w) xi = src.w - 1;
            if (yi >= src.h) yi = src.h - 1;
            for(int c = 0; c < image.channels; ++c){
                PXL_AT(image, x, y, c) = PXL_AT(src, xi, yi, c);
            }
        }
    }

    return image;
}

image_t with_all(double value, int w, int h, int channels){
    image_t i = create_image(w, h, channels);
    FOREACH_PXL(i, {
        PXL_AT(i, x, y, c) = value;
    });
    return i;
}

image_t zeros(int w, int h, int channels){
    return with_all(0, w, h, channels);
}

image_t median(image_t image, int ksize){
    assert(ksize%2 != 0);

    int kend = ksize/2;
    int kstart = -ksize/2;
    
    int middle = (ksize*ksize)/2;
    double_arr array = DOUBLE_ARR(ksize*ksize+1);
    memset(array.items, 0, array.capacity * sizeof(*array.items));
    double aux = 0;
    
    image_t new_image = create_image(image.w, image.h, image.channels);

    int tx = 0, ty = 0;
    FOREACH_PXL(image, {
        for(int kx = kstart; kx <= kend; ++kx){
            for(int ky = kstart; ky <= kend; ++ky){
                tx = x + kx;
                ty = y + ky;
                if(ty >= image.h || ty < 0 || tx < 0 || tx >= image.w){
                    continue;
                }
                APPEND(array, PXL_AT(image, tx, ty, c));
            }
        }
        for(int k = 0; k < array.count; ++k){
            for(int l = 0; l < array.count-1-k; ++l){
                if(VGET(array, l) > VGET(array, l+1)){
                    aux = VGET(array, l);
                    VGET(array, l) = VGET(array, l+1);
                    VGET(array, l+1) = aux;
                }
            }
        }
        PXL_AT(new_image, x, y, c) = VGET(array, middle);
        array.count = 0;
    });
    return new_image;
}

image_t mean(image_t source, int ksize){
    return convolve(source, mean_kernel(ksize));
}

image_t gaussian(image_t source, int ksize, double sigma){
    return convolve(source, gaussian_kernel(ksize, sigma));
}

image_t brightness(image_t image, double b){
    image_t new_image = clone(image);
    FOREACH_PXL(new_image, {
        if(new_image.channels == 4 && c == 3) continue;
        if(PXL_AT(new_image, x, y, c) == VOID_PIXEL) continue;

        PXL_AT(new_image, x, y, c) = b + PXL_AT(new_image, x, y, c);
    });
    return new_image;
}

image_t brightness_channel(image_t image, int channel, double b){
    if(image.channels <= channel) return image;
    for(int y = 0; y < image.h; ++y){
        for(int x = 0; x < image.w; ++x){
            if(PXL_AT(image, x, y, channel) == VOID_PIXEL) continue;
            PXL_AT(image, x, y, channel) = b + PXL_AT(image, x, y, channel);
        }
    }
    return image;
}

image_t brightness_ch(image_t image, int channels, ...){
    if(image.channels < channels) return image;

    va_list args;
    va_start(args, channels);

    double_arr ch_brightness = DOUBLE_ARR(channels);
    for(int c = 0; c < channels; ++c){
        double value = va_arg(args, double);
        APPEND(ch_brightness, value);
    }
    
    for(int y = 0; y < image.h; ++y){
        for(int x = 0; x < image.w; ++x){
            for(int c = 0; c < channels; ++c){
                if(PXL_AT(image, x, y, c) == VOID_PIXEL) continue;
                PXL_AT(image, x, y, c) = ch_brightness.items[c] + PXL_AT(image, x, y, c);
            }
        }
    }
    va_end(args);
    return image;
}

double img_min(image_t image){
    double min = LONG_MAX;
    for(int i = 0; i < image.h*image.w*image.channels; ++i)
        if(image.data[i] < min) min = image.data[i];
    return min;
}

double img_max(image_t image){
    double max = 0;
    for(int i = 0; i < image.h*image.w*image.channels; ++i)
        if(image.data[i] > max) max = image.data[i];
    return max;
}

image_t negative(image_t image){
    return TRANSFORM(clone(image), 255 - pixel);
}

image_t tlog(image_t image, double c){
    c = (c) * (255.0 / log(1.0 + 255.0));
    return TRANSFORM(image, c * log(1.0 + pixel));
}

image_t tpower(image_t image, double gamma, double c){
    return TRANSFORM(image, c * pow(pixel/255.0, gamma) * 255.0);
}

image_t contrast_stretch(image_t image){
    double min = img_min(image);
    double max = img_max(image);
    double diff = max - min;
    return TRANSFORM(image, ((pixel - min)*255) / diff);
}

image_t create_greater_image_between(image_t im1, image_t im2){
    int max_w, min_w, min_h, max_h, max_ch;
    max_w = im1.w;
    max_h = im1.h;
    min_w = im1.w;
    min_h = im1.h;
    max_ch = im1.channels;

    if(im2.channels > max_ch) max_ch = im2.channels;

    if(im2.w > max_w) max_w = im2.w;
    if(im2.h > max_h) max_h = im2.h;

    if(im2.w < min_w) min_w = im2.w;
    if(im2.h < min_h) min_h = im2.h;

    if(im2.h < min_h) min_h = im2.h;

    return create_image(max_w, max_h, max_ch);
}

image_t sum_offset_a(image_t im1, image_t im2, int off_x, int off_y, int ignore_alpha){
    image_t new_image = create_greater_image_between(im1, im2);
    double value = 0;
    int x2, y2;
    FOREACH_PXL(new_image, {
        value = 0;
        x2 = x - off_x;
        y2 = y - off_y;
        if(ignore_alpha && c == 3) PXL_AT(new_image, x, y, c) = 255;
        if(c < im1.channels && (x < im1.w && y < im1.h)) {
            if(PXL_AT(im1, x, y, c) != VOID_PIXEL)
                value += PXL_AT(im1, x, y, c);
        }
        if((x2 >= 0) && (y2 >= 0) && c < im2.channels && (x2 < im2.w && y2 < im2.h)){
            if(PXL_AT(im2, x2, y2, c) != VOID_PIXEL)
                value += PXL_AT(im2, x2, y2, c);
        }
        PXL_AT(new_image, x, y, c) = value;
    });
    return new_image;
}

image_t sum_offset(image_t im1, image_t im2, int off_x, int off_y){
    return sum_offset_a(im1, im2, off_x, off_y, 1);
}

image_t sum(image_t im1, image_t im2){
    return sum_offset(im1, im2, 0, 0);
}

image_t sub_offset_a(image_t im1, image_t im2, int off_x, int off_y, int ignore_alpha){
    image_t new_image = create_greater_image_between(im1, im2);
    double value = 0;
    int x2, y2;
    FOREACH_PXL(new_image, {
        value = 0;
        x2 = x - off_x;
        y2 = y - off_y;
        if(ignore_alpha && c == 3) PXL_AT(new_image, x, y, c) = 255;
        if(c < im1.channels && (x < im1.w && y < im1.h)) {
            if(PXL_AT(im1, x, y, c) != VOID_PIXEL)
                value += PXL_AT(im1, x, y, c);
        }
        if((x2 >= 0) && (y2 >= 0) && c < im2.channels && (x2 < im2.w && y2 < im2.h)){
            if(PXL_AT(im2, x2, y2, c) != VOID_PIXEL)
                value -= PXL_AT(im2, x2, y2, c);
        }
        PXL_AT(new_image, x, y, c) = value;
    });
    return new_image;
}

image_t sub_offset(image_t im1, image_t im2, int off_x, int off_y){
    return sub_offset_a(im1, im2, off_x, off_y, 1);
}

image_t sub(image_t im1, image_t im2){
    return sub_offset(im1, im2, 0, 0);
}

image_t multiply_offset_a(image_t im1, image_t im2, int off_x, int off_y, int ignore_alpha){
    image_t new_image = create_greater_image_between(im1, im2);
    double value = 0;
    int x2, y2;
    FOREACH_PXL(new_image, {
        value = 0;
        x2 = x - off_x;
        y2 = y - off_y;
        if(ignore_alpha && c == 3) PXL_AT(new_image, x, y, c) = 255.0;
        if(c < im1.channels && (x < im1.w && y < im1.h)) {
            if(PXL_AT(im1, x, y, c) != VOID_PIXEL)
                value += PXL_AT(im1, x, y, c);
        }
        if((x2 >= 0) && (y2 >= 0) && c < im2.channels && (x2 < im2.w && y2 < im2.h)){
            if(PXL_AT(im2, x2, y2, c) != VOID_PIXEL)
                value = ((value) * (PXL_AT(im2, x2, y2, c)))/255.0;
        }
        PXL_AT(new_image, x, y, c) = value;
    });
    return new_image;
}

image_t multiply_offset(image_t im1, image_t im2, int off_x, int off_y){
    return multiply_offset_a(im1, im2, off_x, off_y, 1);
}

image_t multiply(image_t im1, image_t im2){
    return multiply_offset(im1, im2, 0, 0);
}

image_t multiply_self(image_t src, int times){
    if(times <= 1) return src;
    image_t result = multiply(src, src);
    image_t aux;
    times -= 1;
    for(int i = 1; i < times; ++i){
        aux = multiply(result, result);
        free_image(&result);
        result = aux;
    }
    return result;
}

image_t blend_offset_a(image_t im1, image_t im2, double b_alpha, int off_x, int off_y, int ignore_alpha){
    image_t new_image = create_greater_image_between(im1, im2);
    double v1 = 0, v2 = 0;
    int x2, y2;
    double alpha = 1 - b_alpha;
    FOREACH_PXL(new_image, {
        x2 = x - off_x;
        y2 = y - off_y;
        v1 = -1;
        v2 = -1;
        if(ignore_alpha && c == 3) PXL_AT(new_image, x, y, c) = 255;
        if(c < im1.channels && (x < im1.w && y < im1.h)) {
            if(PXL_AT(im1, x, y, c) != VOID_PIXEL)
                v1 = PXL_AT(im1, x, y, c);
        }
        if((x2 >= 0) && (y2 >= 0) && c < im2.channels && (x2 < im2.w && y2 < im2.h)){
            if(PXL_AT(im2, x2, y2, c) != VOID_PIXEL)
                v2 = PXL_AT(im2, x2, y2, c);
        }
        if(v1 == -1) v1 = v2;
        if(v2 == -1) v2 = v1;
        if(v2 == -1) {
            v1 = 0;
            v2 = 0;
        }
        PXL_AT(new_image, x, y, c) = (alpha * v1) + (b_alpha * v2);
    });
    return new_image;
}

image_t blend_offset(image_t im1, image_t im2, double b_alpha, int off_x, int off_y){
    return blend_offset_a(im1, im2, b_alpha, off_x, off_y, 1);
}

image_t blend(image_t im1, image_t im2, double b_alpha){
    return blend_offset(im1, im2, b_alpha, 0, 0);
}

image_t screen_offset_a(image_t im1, image_t im2, int off_x, int off_y, int ignore_alpha){
    image_t new_image = create_greater_image_between(im1, im2);
    double v1 = 0, v2 = 0;
    int x2, y2;
    FOREACH_PXL(new_image, {
        x2 = x - off_x;
        y2 = y - off_y;
        v1 = -1;
        v2 = -1;
        if(ignore_alpha && c == 3) PXL_AT(new_image, x, y, c) = 255;
        if(c < im1.channels && (x < im1.w && y < im1.h)) {
            if(PXL_AT(im1, x, y, c) != VOID_PIXEL)
                v1 = PXL_AT(im1, x, y, c);
        }
        if((x2 >= 0) && (y2 >= 0) && c < im2.channels && (x2 < im2.w && y2 < im2.h)){
            if(PXL_AT(im2, x2, y2, c) != VOID_PIXEL)
                v2 = PXL_AT(im2, x2, y2, c);
        }
        if(v1 == -1 && v2 != -1) PXL_AT(new_image, x, y, c) = v2;
        else if(v1 == -1) PXL_AT(new_image, x, y, c) = 0;
        else PXL_AT(new_image, x, y, c) = 255.0 - ((255.0 - v1)*(255.0 - v2))/255;
    });
    return new_image;
}

image_t screen_offset(image_t im1, image_t im2, int off_x, int off_y){
    return screen_offset_a(im1, im2, off_x, off_y, 1);
}
image_t screen(image_t im1, image_t im2){
    return screen_offset(im1, im2, 0, 0);
}

image_t mask_offset_a(image_t src, image_t mask_src, int off_x, int off_y, int ignore_alpha){
    image_t new_image = clone(src);
    int inside_mask;
    int x2, y2;
    FOREACH_PXL(new_image, {
        x2 = x - off_x;
        y2 = y - off_y;
        if(ignore_alpha && c == 3) PXL_AT(new_image, x, y, c) = 255;
        inside_mask = (x2 >= 0) && (y2 >= 0) && c < mask_src.channels && (x2 < mask_src.w && y2 < mask_src.h);

        if(inside_mask && PXL_AT(mask_src, x2, y2, c) == VOID_PIXEL)
            inside_mask = 0;

        if(!inside_mask) PXL_AT(new_image, x, y, c) = 0;
    });
    return new_image;
}

image_t mask_offset(image_t src, image_t mask_src, int off_x, int off_y){
    return mask_offset_a(src, mask_src, off_x, off_y, 1);
}
image_t mask(image_t src, image_t mask_src){
    return mask_offset(src, mask_src, 0, 0);
}

image_t and_offset_a(image_t im1, image_t im2, int off_x, int off_y, int ignore_alpha){
    image_t new_image = create_greater_image_between(im1, im2);
    unsigned char v1 = 0, v2 = 0;
    int x2, y2;
    FOREACH_PXL(new_image, {
        x2 = x - off_x;
        y2 = y - off_y;
        if(ignore_alpha && c == 3) PXL_AT(new_image, x, y, c) = 255;

        if(c < im1.channels && (x < im1.w && y < im1.h)){
            v1 = CLAMP(PXL_AT(im1, x, y, c), 0, 255);
        }
        else v1 = 0;
        if((x2 >= 0) && (y2 >= 0) && c < im2.channels && (x2 < im2.w && y2 < im2.h)){
            v2 = CLAMP(PXL_AT(im2, x2, y2, c), 0, 255);
        }
        else v2 = 0;
        PXL_AT(new_image, x, y, c) = v1 & v2;
    });
    return new_image;
}

image_t and_offset(image_t im1, image_t im2, int off_x, int off_y){
    return and_offset_a(im1, im2, off_x, off_y, 1);
}
image_t and(image_t im1, image_t im2){
    return and_offset(im1, im2, 0, 0);
}


image_t xor_offset_a(image_t im1, image_t im2, int off_x, int off_y, int ignore_alpha){
    image_t new_image = create_greater_image_between(im1, im2);
    unsigned char v1 = 0, v2 = 0;
    int x2, y2;
    FOREACH_PXL(new_image, {
        x2 = x - off_x;
        y2 = y - off_y;
        if(ignore_alpha && c == 3) PXL_AT(new_image, x, y, c) = 255;

        if(c < im1.channels && (x < im1.w && y < im1.h)){
            v1 = CLAMP(PXL_AT(im1, x, y, c), 0, 255);
        }
        else v1 = 0;
        if((x2 >= 0) && (y2 >= 0) && c < im2.channels && (x2 < im2.w && y2 < im2.h)){
            v2 = CLAMP(PXL_AT(im2, x2, y2, c), 0, 255);
        }
        else v2 = 0;
        PXL_AT(new_image, x, y, c) = v1 ^ v2;
    });
    return new_image;
}

image_t xor_offset(image_t im1, image_t im2, int off_x, int off_y){
    return xor_offset_a(im1, im2, off_x, off_y, 1);
}
image_t xor(image_t im1, image_t im2){
    return xor_offset(im1, im2, 0, 0);
}

image_t fill_all(image_t src, double value){
    image_t new_image = clone(src);
    FOREACH_PXL(new_image, {
        if(PXL_AT(new_image, x, y, c) != VOID_PIXEL)
            PXL_AT(new_image, x, y, c) = value;
    });
    return new_image;
}

image_t fill(image_t src, int x, int y, double value){
    image_t new_image = clone(src);
    char * visited = ALLOC(new_image.w * new_image.h * sizeof(*visited));
    memset(visited, 0, new_image.w * new_image.h);

    vec2_arr positions = VEC2_ARR(10);
    PUSH(positions, vec2(x, y));
    
    int match = 0;
    int dx[] = {0, -1, 1, 0};
    int dy[] = {-1, 0, 0, 1};
    vec2_t curr;
    int i = 0, nx, ny;
    
    while(positions.count > 0){
        DEQUEUE(positions, curr);
        if(curr.x >= new_image.w || curr.y >= new_image.h || curr.x < 0 || curr.y < 0)
            continue;

        i = curr.y * new_image.w + curr.x;
        if(visited[i]) continue;
        visited[i] = 1;

        match = 1;
        for(int c = 0; c < new_image.channels; c++){
            if(PXL_AT(src, curr.x, curr.y, c) != PXL_AT(src, x, y, c)){
                match = 0;
                break;
            }
        }
        if(!match) continue;

        for(int c = 0; c < new_image.channels; c++)
            PXL_AT(new_image, curr.x, curr.y, c) = value;
    
        for(int i = 0; i < 4; i++) {
            nx = curr.x + dx[i];
            ny = curr.y + dy[i];
            PUSH(positions, vec2(nx, ny));
        }
    }

    FREE(visited);
    return new_image;
}


image_t flood_fill(image_t src, int x, int y, hex_t color){
    image_t new_image = clone(src);
    char * visited = ALLOC(new_image.w * new_image.h * sizeof(*visited));
    memset(visited, 0, new_image.w * new_image.h);

    vec2_arr positions = VEC2_ARR(10);
    PUSH(positions, vec2(x, y));
    
    int match = 0;
    int dx[] = {0, -1, 1, 0};
    int dy[] = {-1, 0, 0, 1};
    vec2_t curr;
    int i = 0, nx, ny;

    int a = (color >> 24) & 0xFF;
    int r = (color >> 16) & 0xFF;
    int g = (color >> 8) & 0xFF;
    int b = color & 0xFF;
    color = (a << 24) | (b << 16) | (g << 8) | (r);

    while(positions.count > 0){
        DEQUEUE(positions, curr);
        if(curr.x >= new_image.w || curr.y >= new_image.h || curr.x < 0 || curr.y < 0)
            continue;

        i = curr.y * new_image.w + curr.x;
        if(visited[i]) continue;
        visited[i] = 1;

        match = 1;
        for(int c = 0; c < new_image.channels; c++){
            if(PXL_AT(src, curr.x, curr.y, c) != PXL_AT(src, x, y, c)){
                match = 0;
                break;
            }
        }
        if(!match) continue;

        for(int c = 0; c < new_image.channels; c++){
            PXL_AT(new_image, curr.x, curr.y, c) = (0xFF & (color >> (c * 8)));
        }
    
        for(int i = 0; i < 4; i++) {
            nx = curr.x + dx[i];
            ny = curr.y + dy[i];
            PUSH(positions, vec2(nx, ny));
        }
    }

    FREE(visited);
    return new_image;
}

void histdump_from_data(int * histogram) {

    int max_freq = 0;
    for(int i = 0; i < 256; ++i) {
        if(histogram[i] > max_freq)
            max_freq = histogram[i];
    }

    FILE *gnuplot = popen("gnuplot -persistent", "w");
    if (gnuplot == NULL) {
        perror("Failed to open gnuplot");
        exit(EXIT_FAILURE);
    }

    fprintf(gnuplot, "set terminal qt size 800,500 enhanced font 'Verdana,10'\n");
    fprintf(gnuplot, "set title 'Image Histogram' font ',14' textcolor '#333333'\n");
    fprintf(gnuplot, "set xlabel 'Pixel Intensity (0-255)' font ',12' textcolor '#333333'\n");
    fprintf(gnuplot, "set ylabel 'Frequency' font ',12' textcolor '#333333'\n");
    fprintf(gnuplot, "set xrange [0:255]\n");
    // fprintf(gnuplot, "set yrange [0:*]\n");
    // fprintf(gnuplot, "set yrange [0:%.0f]\n", max_freq+10);
    
    fprintf(gnuplot, "set style fill solid 0.8 border rgb '#333333'\n");
    fprintf(gnuplot, "set boxwidth 0.9 relative\n");
    fprintf(gnuplot, "set grid xtics ytics lc rgb '#dddddd' lw 1\n");
    fprintf(gnuplot, "set tics font ',10'\n");
    fprintf(gnuplot, "set border 3 back lc rgb '#666666' lw 2\n");
    fprintf(gnuplot, "set key off\n");
    fprintf(gnuplot, "set object 1 rectangle from screen 0,0 to screen 1,1 behind fillcolor rgb '#f8f8f8' fillstyle solid\n");
    
    fprintf(gnuplot, "set palette defined (0 '#5e81ac', 1 '#81a1c1', 2 '#88c0d0', 3 '#8fbcbb')\n");
    
    fprintf(gnuplot, "plot '-' using 1:($2):1 with boxes lc palette\n");

    for (int i = 0; i < 256; i++) {
        fprintf(gnuplot, "%d %d\n", i, histogram[i]);
    }
    fprintf(gnuplot, "e\n");
    fflush(gnuplot);
    pclose(gnuplot);
}

void histdump(image_t image) {
    double_arr frequencies = DOUBLE_ARR(256);
    for(int i = 0; i < frequencies.capacity; ++i) 
        frequencies.items[i] = 0.0;

    double value = 0.0;
    int idx = 0;
    for(int i = 0; i < image.h; ++i) {
        for(int j = 0; j < image.w; ++j) {
            value = PXL_AT(image, j, i, 0);
            idx = (int)value;
            idx = (idx > 255) ? 255 : ((idx < 0) ? 0 : idx);
            frequencies.items[idx] += 1.0;
        }
    }

    double max_freq = 0.0;
    for(int i = 0; i < 256; ++i) {
        if(frequencies.items[i] > max_freq)
            max_freq = frequencies.items[i];
    }

    FILE *gnuplot = popen("gnuplot -persistent", "w");
    if (gnuplot == NULL) {
        perror("Failed to open gnuplot");
        exit(EXIT_FAILURE);
    }

    fprintf(gnuplot, "set terminal qt size 800,500 enhanced font 'Verdana,10'\n");
    fprintf(gnuplot, "set title 'Image Histogram' font ',14' textcolor '#333333'\n");
    fprintf(gnuplot, "set xlabel 'Pixel Intensity (0-255)' font ',12' textcolor '#333333'\n");
    fprintf(gnuplot, "set ylabel 'Frequency' font ',12' textcolor '#333333'\n");
    fprintf(gnuplot, "set xrange [0:255]\n");
    // fprintf(gnuplot, "set yrange [0:90]\n");
    // fprintf(gnuplot, "set yrange [0:%.0f]\n", max_freq+10);
    
    fprintf(gnuplot, "set style fill solid 0.8 border rgb '#333333'\n");
    fprintf(gnuplot, "set boxwidth 0.9 relative\n");
    fprintf(gnuplot, "set grid xtics ytics lc rgb '#dddddd' lw 1\n");
    fprintf(gnuplot, "set tics font ',10'\n");
    fprintf(gnuplot, "set border 3 back lc rgb '#666666' lw 2\n");
    fprintf(gnuplot, "set key off\n");
    fprintf(gnuplot, "set object 1 rectangle from screen 0,0 to screen 1,1 behind fillcolor rgb '#f8f8f8' fillstyle solid\n");
    
    fprintf(gnuplot, "set palette defined (0 '#5e81ac', 1 '#81a1c1', 2 '#88c0d0', 3 '#8fbcbb')\n");
    
    fprintf(gnuplot, "plot '-' using 1:($2):1 with boxes lc palette\n");

    for (int i = 0; i < 256; i++) {
        fprintf(gnuplot, "%d %.2f\n", i, frequencies.items[i]);
    }
    fprintf(gnuplot, "e\n");
    fflush(gnuplot);
    pclose(gnuplot);
}