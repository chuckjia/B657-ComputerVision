//
// detect.cpp : Detect integrated circuits in printed circuit board (PCB) images.
//
// Based on skeleton code by D. Crandall, Spring 2018
//
// PUT YOUR NAMES HERE
//
//

#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent squared gradient magnitudes,
// harris corner scores, etc. 
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in 
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below is a helper functions that overlays rectangles
// on an image plane for visualization purpose. 

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//
void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
	for(int w=-width/2; w<=width/2; w++) {
		int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;

		// if any of the coordinates are out-of-bounds, truncate them
		top = min( max( top, 0 ), input.rows()-1);
		bottom = min( max( bottom, 0 ), input.rows()-1);
		left = min( max( left, 0 ), input.cols()-1);
		right = min( max( right, 0 ), input.cols()-1);

		// draw top and bottom lines
		for(int j=left; j<=right; j++)
			input[top][j] = input[bottom][j] = graylevel;
		// draw left and right lines
		for(int i=top; i<=bottom; i++)
			input[i][left] = input[i][right] = graylevel;
	}
}

// DetectedBox class may be helpful!
//  Feel free to modify.
//
class DetectedBox {
public:
	int row, col, width, height;
	double confidence;
};

// Function that outputs the ascii detection output file
void  write_detection_txt(const string &filename, const vector<DetectedBox> &ics)
{
	ofstream ofs(filename.c_str());

	for(vector<DetectedBox>::const_iterator s=ics.begin(); s != ics.end(); ++s)
		ofs << s->row << " " << s->col << " " << s->width << " " << s->height << " " << s->confidence << endl;
}

// Function that outputs a visualization of detected boxes
void  write_detection_image(const string &filename, const vector<DetectedBox> &ics, const SDoublePlane &input)
{
	SDoublePlane output_planes[3];

	for(int p=0; p<3; p++)
	{
		output_planes[p] = input;
		for(vector<DetectedBox>::const_iterator s=ics.begin(); s != ics.end(); ++s)
			overlay_rectangle(output_planes[p], s->row, s->col, s->row+s->height-1, s->col+s->width-1, p==2?255:0, 2);
	}

	SImageIO::write_png_file(filename.c_str(), output_planes[0], output_planes[1], output_planes[2]);
}



// The rest of these functions are incomplete. These are just suggestions to 
// get you started -- feel free to add extra functions, change function
// parameters, etc.

// Convolve an image with a separable convolution kernel
//
SDoublePlane convolve_separable(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{
	SDoublePlane output(input.rows(), input.cols());

	// Convolution code here

	return output;
}

// Convolve an image with a  convolution kernel
//
SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter)
{
	SDoublePlane output(input.rows(), input.cols());

	// Convolution code here
	int nrow_im = input.rows(), ncol_im = input.cols(), nrow_f = filter.rows(), ncol_f = filter.cols();

	return output;
}


// Apply a sobel operator to an image, returns the result
// 
SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx)
{
	SDoublePlane output(input.rows(), input.cols());

	// Implement a sobel gradient estimation filter with 1-d filters


	return output;
}

// Apply an edge detector to an image, returns the binary edge map
// 
SDoublePlane find_edges(const SDoublePlane &input, double thresh=0)
{
	SDoublePlane output(input.rows(), input.cols());

	// Implement an edge detector of your choice, e.g.
	// use your sobel gradient operator to compute the gradient magnitude and threshold

	return output;
}

// Pre-process image to eliminate noises
SDoublePlane preprocess(const SDoublePlane &input_r, SDoublePlane &input_g, SDoublePlane &input_b) {
	double lower = 20, upper = 90;
	int nrow = input_r.rows(), ncol = input_r.cols();
	SDoublePlane output(nrow, ncol);
	for (int i = 0; i < nrow; ++i)
		for (int j = 0; j < ncol; ++j) {
			double r = input_r[i][j], g = input_g[i][j], b = input_b[i][j];
			if (r > lower & r < upper & g > lower & g < upper & b > lower & b < upper) {
				output[i][j] = 0;
				// printf("%f\n", output[i][j]);
			} else
				output[i][j] = 255;
		}
	SDoublePlane newoutput = output;
	for (int i = 0; i < nrow; ++i)
		for (int j = 0; j < ncol; ++j)
			if (output[i][j] < 200) {
				int white_count = 0;
				for (int ii = 0; ii < 5; ++ii)
					for (int jj = 0; jj < 5; ++jj)
						white_count += output[i][j] > 200;
				if (white_count > 13)
					newoutput[i][j] = 255;
			}
	return newoutput;
}


//
// This main file just outputs a few test images. You'll want to change it to do
//  something more interesting!
//
int main(int argc, char *argv[])
{
	if(!(argc == 2))
	{
		cerr << "usage: " << argv[0] << " input_image" << endl;
		return 1;
	}

	string input_filename(argv[1]);
	SDoublePlane input_image= SImageIO::read_png_file(input_filename.c_str());
	SDoublePlane input_r(input_image.rows(), input_image.cols()),
			input_g(input_image.rows(), input_image.cols()),
			input_b(input_image.rows(), input_image.cols());
	SImageIO::read_png_file(input_filename.c_str(), input_r, input_g, input_b);
	SDoublePlane output_prep = preprocess(input_r, input_g, input_b);
	string test_filename = "abc.png";
	SImageIO::write_png_file(test_filename.c_str(), output_prep, output_prep, output_prep);

	// test step 2 by applying mean filters to the input image
	SDoublePlane mean_filter(3,3);
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			mean_filter[i][j] = 1/9.0;
	SDoublePlane output_image = convolve_general(input_image, mean_filter);


	// randomly generate some detected ics -- you'll want to replace this
	//  with your ic detection code obviously!
	vector<DetectedBox> ics;
	for(int i=0; i<10; i++)
	{
		DetectedBox s;
		s.row = rand() % input_image.rows();
		s.col = rand() % input_image.cols();
		s.width = 20;
		s.height = 20;
		s.confidence = rand();
		ics.push_back(s);
	}

	write_detection_txt("detected.txt", ics);
	write_detection_image("detected.png", ics, input_image);
}
