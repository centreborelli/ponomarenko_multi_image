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

More options are displayed

## Test
Run the extended version on multiple images:
```
./ponomarenko_multi ../test_images/seq0/
```

Run the original version on a single image:
```
./ponomarenko_single ../test_images/seq0/image-001.png
```

Since the `test_images/seq0` contains only one image, the outputs from two programs should be the same.

Example:

```
yanhao@li:~/ponomarenko_multi_image/build$ ./ponomarenko_multi ../test_images/seq0/
97.556854  20.464422  65.437210  2.571375  3.398506  3.360783  
161.559494  65.373238  77.609665  2.299690  3.288906  3.380521  
176.508362  92.205841  95.019363  2.108210  2.986117  3.977546  
208.514084  108.968460  110.455032  1.544361  2.701997  4.355613  
221.660721  134.267456  125.007042  1.287640  2.412065  4.302766  
240.984161  203.753891  181.675766  0.807139  1.327316  3.156296
```

```
yanhao@li:~/ponomarenko_multi_image/build$ ./ponomarenko_single ../test_images/seq0/image-001.png
97.556854  20.464422  65.437210  2.571375  3.398506  3.360783  
161.559494  65.373238  77.609665  2.299690  3.288906  3.380521  
176.508362  92.205841  95.019363  2.108210  2.986117  3.977546  
208.514084  108.968460  110.455032  1.544361  2.701997  4.355613  
221.660721  134.267456  125.007042  1.287640  2.412065  4.302766  
240.984161  203.753891  181.675766  0.807139  1.327316  3.156296
```