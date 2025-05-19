#include "imagec.h"

int main(){
    MM_BEGIN();

    /**
     * Original Figure & Gaussian Filter
     */
    image_t fig1041 = load_image("./sample/fig1041.png");

    /**
     * Apply sobel operator to get gradient magnitude
     */
    image_t gradient = sobel(fig1041);

    /**
     * Threshold gradient on the 99.7 percentile
     */
    image_t limiarized = threshold_percentile(gradient, 0.997);

    /**
     * Multiply by limiarized
     */
    image_t multiplied = multiply(limiarized, fig1041);

    /**
     * Calculate histogram from product result
     */
    int * histogram_multiplied = histdata__ignore_zeros(multiplied);

    /**
     * Apply Otsu-Threshold with the multiplied histogram calculated
     */
    image_t otsu = otsu_from_hist(fig1041, histogram_multiplied);

    write_image("./Output 1 - original.png", fig1041);  
    write_image("./Output 2 - gradient.png", gradient);
    write_image("./Output 3 - gradient-limiarized.png", limiarized);
    write_image("./Output 4 - multiplied.png", multiplied);
    write_image("./Output 5 - otsu-threshold.png", otsu);

    MM_CLEAN();
    return 0;
}