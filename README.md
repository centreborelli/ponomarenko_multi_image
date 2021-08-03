# ponomarenko_multi_image
Extension of Ponomarenko Noise Estimation method for a sequence of images.

<<<<<<< HEAD
## Usage
=======
## Usage:
>>>>>>> 618404f5931ce61ee51e048f66a8fd33c8909fcc

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

