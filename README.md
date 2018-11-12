# B657 Computer Vision

This repository holds the code for my projects from the course B657 Computer Vision in Spring 2018.

## Project 1 (HW1A)

In this project, we built a C++ program that adds invisible watermarks to photos. The added watermarks, although cannot be seen by human eyes, can be detected by our program. The algorithm uses Fourier transform to add watermarks to the Fourier space.

To compile the source code, use gcc on `watermarks.cpp`
```
g++ watermarks.cpp -o p2
```
To run the program, use the format
```
p2 [problemID] [inputfilename] [outputfilename]
```
The `[problemID]` corresponds to `1.1`, `1.2`, `1.3`, and `1.4`, which are 4 features of the program. See the `.pdf` file for details.

&nbsp;


## Project 2 (HW1B)

This project is a C++ program that detects integrated circuits on printed circuit boards (PCBs). Multiple techniques were used, including Hough transformation based edge detection,
sliding windows, Fourier transform based image similarity analysis, and connected-component labeling.

To compile, use gcc on `detect.cpp`
```
g++ detect.cpp -o detect
```
To run the program, use
```
detect [inputfilename]
```
The output of the program include images of the PCB with boxes around the detected integrated circuits and a `.txt` file detailing the coordinates of the detected areas.

&nbsp;

## Project 3 (HW2)

In this project, a C++ program was created to generate panoramic photos. The program also has the features of blending images and performing perspective changes on photos. The algorithm used techniques including pyramid representations for blending, RANSAC and SIFT for matching and stitching images, and homography transformation for perspective changes.

To compile, use gcc on the `a2.cpp` file. To run the program used
```
a2 [partID] [inputfile1] [inputfile2] [inputfile3]
```

Here `[partID]` can be the following options
* `part1`: Changing the perspective of an image. It can also add the modified image into a prescribed area of another image. For example, we can change the perspective of a headshot and fit it into a billboard in another photo.

* `part2`: Blend 2 images. The blending area would be in the middle

* `part3`: Image matching. This feature can match two images and output an image marking the matched areas between the two images.

* `part4`: Creating panoramic photos using multiple image input.
