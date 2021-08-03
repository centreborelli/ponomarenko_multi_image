# ponomarenko_multi_image
Extension of Ponomarenko Noise Estimation method for a sequence of images.

## Usage

Compile:

```
mkdir build && cd build
cmake ..
make -j
```

Run the noise estimation program for a directory containing `.png` image files:
```
./ponomarenko_multi <image-dir>
```

## Test
Run the extended version on multiple images
```
./ponomarenko_multi ../test_images/seq0/
```

Run the original version on a single image:
```
./ponomarenko_single ../test_images/seq0/image-001.png
```

Since the `test_images/seq0` contains only one image, the outputs from two programs should be the same.

