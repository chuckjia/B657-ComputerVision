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
#include <list>

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
	// Convolution code here
	int nrow = input.rows(), ncol = input.cols();
	int row_ftr_len = row_filter.cols(), col_ftr_len = col_filter.rows();  // Size of the filters

	// Row filter
	SDoublePlane input_interm(nrow, ncol);
	int fcenter = col_ftr_len / 2;
	for (int i = 0; i < nrow; ++i) {
		for (int j = 0; j < ncol; ++j) {
			double sum = 0;
			for (int k_ftr = 0; k_ftr < row_ftr_len; ++k_ftr) {
				int k_img = j + fcenter - k_ftr;
				sum += row_filter[0][k_ftr] *
						((k_img < 0 || k_img >= ncol) ? input[i][j] : input[i][k_img]);
			}
			input_interm[i][j] = sum;
		}
	}
	// Column filter
	SDoublePlane output(nrow, ncol);
	fcenter = col_ftr_len / 2;
	for (int i = 0; i < nrow; ++i) {
		for (int j = 0; j < ncol; ++j) {
			double sum = 0;
			for (int k_ftr = 0 ; k_ftr < col_ftr_len; ++k_ftr) {
				int k_img = i + fcenter - k_ftr;
				sum += col_filter[k_ftr][0] *
						((k_img < 0 || k_img >= nrow) ? input[i][j] : input_interm[k_img][j]);
			}
			output[i][j] = sum;
		}
	}
	return output;
}

// Convolve an image with a  convolution kernel
//
SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter)
{
	// Convolution code here
	int nrow = input.rows(), ncol = input.cols(), nrow_ftr = filter.rows(), ncol_ftr = filter.cols();
	SDoublePlane output(nrow, ncol);

	int fcenter_row = nrow_ftr / 2, fcenter_col = ncol_ftr / 2;
	for (int i = 0; i < nrow; ++i)
		for (int j = 0; j < ncol; ++j) {
			double sum = 0;
			for (int ii = 0 ; ii < nrow_ftr; ++ii)
				for (int jj = 0; jj < ncol_ftr; ++jj) {
					int i_curr = i + fcenter_row - ii, j_curr = j + fcenter_col - jj;
					sum += filter[ii][jj] * (
							(i_curr < 0 || i_curr >= nrow || j_curr < 0 || j_curr >= ncol) ? input[i][j] : input[i_curr][j_curr]);
				}
			output[i][j] += sum;
		}
	return output;
}


// Apply a sobel operator to an image, returns the result
// 
SDoublePlane sobel_gradient_filter(const SDoublePlane &input)
{
	// Implement a sobel gradient estimation filter with 1-d filters
	int nrow = input.rows(), ncol = input.cols();
	SDoublePlane output(nrow, ncol);
	SDoublePlane sobel_x(3,3);
	sobel_x[0][0] = -1; sobel_x[0][1] = 0; sobel_x[0][2] = 1;
	sobel_x[1][0] = -2; sobel_x[1][1] = 0; sobel_x[1][2] = 2;
	sobel_x[2][0] = -1; sobel_x[2][1] = 0; sobel_x[2][2] = 1;

	SDoublePlane sobel_y(3,3);
	sobel_y[0][0] = -1; sobel_y[0][1] = -2; sobel_y[0][2] = -1;
	sobel_y[1][0] = 0; sobel_y[1][1] = 0; sobel_y[1][2] = 0;
	sobel_y[2][0] = 1; sobel_y[2][1] = 2; sobel_y[2][2] = 1;

	// Implement a sobel gradient estimation filter with 1-d filters
	SDoublePlane output_1 = convolve_general(input, sobel_x);
	SDoublePlane output_2 = convolve_general(input, sobel_y);
	for (int i = 0; i < nrow; ++i)
		for (int j = 0; j < ncol; ++j)
			output[i][j] = sqrt(pow(output_1[i][j], 2) + pow(output_2[i][j], 2));
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
	int nrow = input_r.rows(), ncol = input_r.cols();
	SDoublePlane output(nrow, ncol);
	for (int i = 0; i < nrow; ++i)
		for (int j = 0; j < ncol; ++j) {
			double r = input_r[i][j], g = input_g[i][j], b = input_b[i][j];
			if (r < 100 || g < 130 || b < 100) {
				double mean = 1. / 3 * (r + g + b),
						var = 1. / 3 * (pow(r - mean, 2) + pow(g - mean, 2) + pow(b - mean, 2));
				if (var < 150)
					output[i][j] = 0;  // Black
				else
					output[i][j] = 255;  // White
			} else
				output[i][j] = 255;
		}
	SDoublePlane newoutput = output;
	int num_layer = 3, max_i = nrow - num_layer - 1, max_j = ncol - num_layer - 1;
	for (int i = num_layer; i < max_i; ++i)
		for (int j = num_layer; j < max_j; ++j)
			if (output[i][j] < 200) {
				int white_count = 0;
				for (int ii = -num_layer; ii <= num_layer; ++ii)
					for (int jj = -num_layer; jj <= num_layer; ++jj)
						white_count += output[i + ii][j + ii] > 200;
				if (white_count > 0.2 * pow(2 * num_layer + 1, 2))
					newoutput[i][j] = 255;
			}

	// Check every column
	for (int j = 0; j < ncol; ++j) {
		int count = 0;
		for (int i = 0; i < nrow; ++i)
			count += newoutput[i][j] < 100;  // Choose black
		if (count < nrow * 0.1)
			for (int i = 0; i < nrow; ++i)
				newoutput[i][j] = 255;
	}
	// Check every row
	for (int i = 0; i < nrow; ++i) {
		int count = 0;
		for (int j = 0; j < ncol; ++j)
			count += newoutput[i][j] < 100;  // Choose black
		if (count < ncol * 0.1)
			for (int j = 0; j < nrow; ++j)
				newoutput[i][j] = 255;
	}
	return newoutput;
}


_DTwoDimArray<bool> hough_transform(const _DTwoDimArray<bool> &input) {
	const int th_rad = 1,   // Radius of theta interveal. Max thickness of HALF of the edges/lines in terms of number of theta units
			th = 2 * th_rad + 1;  // Max thickness of the edges/lines in terms of theta
	double theta_unit = M_PI / 512;  // Increment of theta

	int nrow = input.rows(), ncol = input.cols(), horiz[th][nrow], vert[th][ncol];
	for (int i = 0; i < th; ++i) {
		for (int j = 0; j < nrow; ++j)  horiz[i][j] = 0;
		for (int j = 0; j < ncol; ++j)  vert[i][j] = 0;
	}

	for (int i = 0; i < nrow; ++i)
		for (int j = 0; j < ncol; ++j) {
			int x = j, y = nrow - 1 - i;
			if (input[y][x])
				for (int kk = -th_rad; kk < th_rad; ++kk) {
					// Horizontal line, around normal theta = PI / 2
					double theta = 0.5 * M_PI + kk * theta_unit;
					int rho = (int) (x * cos(theta) + y * sin(theta));
					if (rho >= 0 && rho < nrow)
						horiz[kk + th_rad][rho] += 1;

					// Vertical line, around theta = PI / 2
					theta = kk * theta_unit;
					rho = (int) (x * cos(theta) + y * sin(theta));
					if (rho >= 0 && rho < ncol)
						vert[kk + th_rad][rho] += 1;
				}
		}

	// Find lines with high votes
	int min_ed_horiz = 200, //(4 / 5.) * th * input.cols(),  // Minimum edge length
			min_ed_vert = 200; // (2 / 3.) * th * input.rows();
	int bin_rad = 2,  // Radius of the bin of rho
			bin_size = 2 * bin_rad + 1;

	// Find horizontal lines and store them in a list of int consisting of rho values
	std::list<int> horiz_lines;  // With normal theta = PI / 2
	for (int rho = 0; rho < nrow; ++rho) {
		int count = 0;
		for (int i = 0; i < th; ++i)  // Iterate over theta - eps to theta + eps
			for (int j = -bin_size; j < bin_size; ++j) // Iterate
				count += (rho + j < 0 && rho + j >= nrow) ? 0 : horiz[i][rho + j];
		if (count > min_ed_horiz)
			horiz_lines.push_back(rho);
	}

	// Find vertical lines
	std::list<int> vert_lines;
	for (int rho = 0; rho < ncol; ++rho) {
		int count = 0;
		for (int i = 0; i < th; ++i)
			for (int j = -bin_size; j < bin_size; ++j)
				count += (rho + j < 0 || rho + j > ncol) ? 0 : vert[i][rho + j];
		if (count > min_ed_vert)
			vert_lines.push_back(rho);
	}

	_DTwoDimArray<bool> output(nrow, ncol);
	memset(output.data_ptr(), 0, sizeof(bool) * nrow * ncol);

	for (std::list<int>::iterator it = horiz_lines.begin(); it != horiz_lines.end(); ++it){
		int row = *it;
		for (int j = 0; j < ncol; ++j)
			output[row][j] = true;
	}

	for (std::list<int>::iterator it = vert_lines.begin(); it != vert_lines.end(); ++it) {
		int col = *it;
		for (int j = 0; j < nrow; ++j)
			output[j][col] = true;
	}
	// std::cout << "Vertical lines: " << *it << "\n";

	return output;
}

SDoublePlane bool_to_doubleplane(_DTwoDimArray<bool> &input) {
	int nrow = input.rows(), ncol = input.cols();
	SDoublePlane output(nrow, ncol);
	for (int i = 0; i < nrow; ++i)
		for (int j = 0; j < ncol; ++j)
			output[i][j] = (input[i][j] ? 255 : 0);
	return output;
}

_DTwoDimArray<bool> doubleplane_to_bool_black(SDoublePlane &input) {
	int nrow = input.rows(), ncol = input.cols();
	_DTwoDimArray<bool> output(nrow, ncol);
	for (int i = 0; i < nrow; ++i)
		for (int j = 0; j < ncol; ++j)
			output[i][j] = input[i][j] < 100;  // True if black
	return output;
}

double calc_true_frac(_DTwoDimArray<bool> &input, int i0, int i1, int j0, int j1) {
	int count = 0;
	for (int i = i0; i <= i1; ++i)
		for (int j = j0; j <= j1; ++j)
			count += input[i][j] ? 1 : 0;
	return ((double) count) / ((i1 - i0 + 1) * (j1 - j0 + 1));
}

_DTwoDimArray<bool> slid_win(SDoublePlane &input) {
	_DTwoDimArray<bool> input_bool = doubleplane_to_bool_black(input);
	_DTwoDimArray<bool> ans(input.rows(), input.cols());
	memset(ans.data_ptr(), 0, sizeof(bool) * ans.rows() * ans.cols());

	int win_size = 5, check_size = 10;  // Check size is the size of the sequence to be checked in 1D subsampling
	double tol = 0.7;  // If the fraction of true points within one window passes tol, then the window will be marked true
	int row_start = win_size, row_end = input.rows() - win_size - 1,
			col_start = win_size, col_end = input.cols() - win_size - 1;
	for (int i = row_start; i < row_end; i += win_size)
		for (int j = col_start; j < col_end; j += win_size)
			if (calc_true_frac(input_bool, i + 1, i + 1, j, j + check_size) > tol &&
					calc_true_frac(input_bool, i, i + check_size, j + 1, j + 1) > tol)  // Check row under the start corner and the col to the right the corner to avoid serrated edges
				for (int k = 0; k < win_size; ++k)
					for (int kk = 0; kk < win_size; ++kk)
						ans[i + k][j + kk] = true;
	return ans;
}

vector<DetectedBox> mark_win(_DTwoDimArray<bool> &input, int win_size) {
	int nrow = input.rows(), ncol = input.cols();
	vector<DetectedBox> ans;
	int row_start = win_size, row_end = input.rows() - win_size - 1,
			col_start = win_size, col_end = input.cols() - win_size - 1;
	int grid_nrow = (row_end - row_start) / win_size, grid_ncol = (col_end - col_start) / win_size;
	_DTwoDimArray<bool> grid_unused(grid_nrow + 2, grid_ncol + 2);
	memset(grid_unused.data_ptr(), 1, sizeof(bool) * grid_unused.rows() * grid_unused.cols());

	for (int i = 0; i < grid_nrow; ++i) {
		int row = row_start + i * win_size;
		// printf("row = %d\n", row);

		for (int j = 0; j < grid_ncol; ++j) {
			int col = col_start + j * win_size;
			if (input[row][col] && grid_unused[i][j]) {  // If a window beginning at (row, col) is marked as true
				int score, score_tol = 2;
				// Check rows below
				int box_curr_row = i + 1, row_curr = row + win_size, last_row = row;
				score = 0;
				while (row_curr < row_end && score <= score_tol) {
					if (!input[row_curr][col + win_size] || !grid_unused[box_curr_row][j + 1])  // If a window is false or used
						++score;
					else  // If a window is true
						last_row = row_curr;
					row_curr += win_size;
					++box_curr_row;
				}

				// Check cols to the right
				int box_curr_col = j + 1, col_curr = col + win_size, last_col = col;
				score = 0;
				while (col_curr < col_end && score <= score_tol) {
					if (!input[row + win_size][col_curr] || !grid_unused[i + 1][box_curr_col])
						++score;
					else
						last_col = col_curr;
					col_curr += win_size;
					++box_curr_col;
				}

				DetectedBox box;
				box.row = row - win_size;
				box.col = col - win_size;
				box.height = last_row - row + 2 * win_size;
				box.width = last_col - col + 2 * win_size;
				box.confidence = 1;
				ans.push_back(box);

				int last_box_index_row = (last_row - row_start) / win_size + 1,
						last_box_index_col = (last_col - col_start) / win_size + 1;
				for (int ii = i; ii <= last_box_index_row; ++ii)
					for (int jj = j; jj <= last_box_index_col; ++jj)
						grid_unused[ii][jj] = false;
			}
		}
	}
	return ans;
}


void write_grayscale_png(const char* filename, const SDoublePlane &input) {
	SImageIO::write_png_file(filename, input, input, input);
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

	/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
	 * Test 1: Preprocess
	 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

	string IC_num = input_filename.substr(3, 1);
	string test_filename;
	int nrow = input_image.rows(), ncol = input_image.cols();
	SDoublePlane input_r(nrow, ncol), input_g(nrow, ncol), input_b(nrow, ncol);
	SImageIO::read_png_file(input_filename.c_str(), input_r, input_g, input_b);
	SDoublePlane output_1 = preprocess(input_r, input_g, input_b);

	test_filename = "TT_IC"; test_filename.append(IC_num); test_filename.append("_01prep.png");
	SImageIO::write_png_file(test_filename.c_str(), output_1, output_1, output_1);

	/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
	 * Test 2: Sobel filter
	 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

	/*test_filename = "TT_IC"; test_filename.append(IC_num); test_filename.append("_02sobel.png");
	output_1 = sobel_gradient_filter(output_1);
	write_grayscale_png(test_filename.c_str(), output_1);*/

	/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
	 * Test 3: Hough transform
	 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

	/*_DTwoDimArray<bool> output_bool_1(nrow, ncol);
	for (int i = 0; i < nrow; ++i)
		for (int j = 0; j < ncol; ++j)
			output_bool_1[i][j] = output_1[i][j] > 20;

	output_bool_1 = hough_transform(output_bool_1);
	output_1 = bool_to_doubleplane(output_bool_1);

	test_filename = "TT_IC"; test_filename.append(IC_num); test_filename.append("_03hough.png");
	write_grayscale_png(test_filename.c_str(), output_1);*/

	/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
	 * Test 4: Use middle line to represent beams
	 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

	/* SDoublePlane output2(nrow, ncol);

	// Intercept of horizontal lines / iterate along column 0
	int first = 0;
	bool accum = false;
	for (int i = 0; i < nrow; ++i)
		if (output_1[i][0] > 200)  // White
			if (!accum) {
				first = i;
				accum = true;
			} else { }
		else { // Black
			if (accum) {
				// output2[(first + i - 1) / 2][0] = 255;
				output2[first][0] = 255;
				output2[i - 1][0] = 255;
				first = i;
				accum = false;
			} else
				++first;
		}

	// Intercept of vertical lines / iterate along row 0
	first = 0;
	accum = false;
	for (int i = 0; i < ncol; ++i)
		if (output_1[0][i] > 100)  // White
			if (!accum) {
				first = i;
				accum = true;
			} else {}
		else {  // Black
			if (accum) {
				// output2[0][(first + i - 1) / 2] = 255;
				output2[0][first] = 255;
				output2[0][i - 1] = 255;
				first = i;
				accum = false;
			} else
				++first;
		}

	for (int i = 0; i < nrow; ++i)
		if (output2[i][0] > 200)
			for (int j = 0; j < ncol; ++j)
				output2[i][j] = 255;

	for (int j = 0; j < ncol; ++j)
		if (output2[0][j] > 200)
			for (int i = 0; i < nrow; ++i)
				output2[i][j] = 255;

	test_filename = "TT_IC"; test_filename.append(IC_num); test_filename.append("_04shrunkLines.png");
	write_grayscale_png(test_filename.c_str(), output2); */

	/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
	 * Test 5
	 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

	_DTwoDimArray<bool> output_2_bool = slid_win(output_1);
	SDoublePlane output_slid_win = bool_to_doubleplane(output_2_bool);
	test_filename = "TT_IC"; test_filename.append(IC_num); test_filename.append("_05slidwin.png");
	write_grayscale_png(test_filename.c_str(), output_slid_win);

	// test step 2 by applying mean filters to the input image
	/*SDoublePlane mean_filter(3,3);
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			mean_filter[i][j] = 1/9.0;
	SDoublePlane output_image = convolve_general(input_image, mean_filter);*/


	// randomly generate some detected ics -- you'll want to replace this
	//  with your ic detection code obviously!
	vector<DetectedBox> ics = mark_win(output_2_bool, 5);
	/*for(int i=0; i<10; i++)
	{
		DetectedBox s;
		s.row = rand() % input_image.rows();
		s.col = rand() % input_image.cols();
		s.width = 20;
		s.height = 20;
		s.confidence = rand();
		ics.push_back(s);
	}*/

	write_detection_txt("detected.txt", ics);
	write_detection_image("detected.png", ics, input_image);
}
