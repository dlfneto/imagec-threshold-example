### Using Image Segmentation Techniques In C

```c
#include "imagec.h"

int main(){
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

    return 0;
}
```

Output

```c
write_image("./Output 1 - original.png", fig1041);  
```

<div align="center">
  <img height="270" src="./Output 1 - original.png" /><br/><br/>
</div>

```c
write_image("./Output 2 - gradient.png", gradient);
```
<div align="center">
  <img height="270" src="./Output 2 - gradient.png" /><br/><br/>
</div>

```c
write_image("./Output 3 - gradient-limiarized.png", limiarized);
```
<div align="center">
  <img height="270" src="./Output 3 - gradient-limiarized.png" /><br/><br/>
</div>

```c
write_image("./Output 4 - multiplied.png", multiplied);
```

<div align="center">
  <img height="270" src="./Output 4 - multiplied.png" /><br/><br/>
</div>

```c
write_image("./Output 5 - otsu-threshold.png", otsu);
```
<div align="center">
  <img height="270" src="./Output 5 - otsu-threshold.png" /><br/><br/>
</div>
