//
// Watermark.cpp : Add watermark to an image, or inspect if a watermark is present.
//
// Based on skeleton code by D. Crandall, Spring 2018
//
// PUT YOUR NAMES HERE
//
//

//Link to the header file
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <SImage.h>
#include <SImageIO.h>
#include <fft.h>

using namespace std;

const bool normalize_intensity = false;
int l_CONST = 32;  // Number of bins
double alpha_CONST = 45;
double radius_fraction_CONST = 0.9;

// This code requires that input be a *square* image, and that each dimension
//  is a power of 2; i.e. that input.width() == input.height() == 2^k, where k
//  is an integer. You'll need to pad out your image (with 0's) first if it's
//  not a square image to begin with. (Padding with 0's has no effect on the FT!)
//
// Forward FFT transform: take input image, and return real and imaginary parts.
//
void fft(const SDoublePlane &input, SDoublePlane &fft_real, SDoublePlane &fft_imag)
{
	fft_real = input;
	fft_imag = SDoublePlane(input.rows(), input.cols());

	FFT_2D(1, fft_real, fft_imag);
}

// Inverse FFT transform: take real and imaginary parts of fourier transform, and return
//  real-valued image.
//
void ifft(const SDoublePlane &input_real, const SDoublePlane &input_imag, SDoublePlane &output_real)
{
	output_real = input_real;
	SDoublePlane output_imag = input_imag;

	FFT_2D(0, output_real, output_imag);
}

// Calculate the log of the norm of a complex number real + i * imag
double calc_log_magnitude(double real, double imag) {
	return log(sqrt(real * real + imag * imag));
}

// Write this in Part 1.1
SDoublePlane fft_magnitude(const SDoublePlane &fft_real, const SDoublePlane &fft_imag)
{
	int nrow = fft_real.rows(), ncol = fft_real.cols();
	SDoublePlane fft_mag(nrow, ncol);
	for (int i = 0; i < nrow; i++)
		for (int j = 0; j < ncol; j++)
			fft_mag[i][j] = calc_log_magnitude(fft_real[i][j], fft_imag[i][j]);
	return fft_mag;
}

double get_max_intensity(const SDoublePlane &input) {
	int nrow = input.rows(), ncol = input.cols();
	double max_intensity = -100;
	for (int i = 0; i < nrow; ++i)
		for (int j = 0; j < ncol; ++j)
			if (input[i][j] > max_intensity)
				max_intensity = input[i][j];
	return max_intensity;
}

SDoublePlane normalize_spec_img(const SDoublePlane &input) {
	int nrow = input.rows(), ncol = input.cols();
	double max_intensity = get_max_intensity(input);
	SDoublePlane output(nrow, ncol);
	for (int i = 0; i < nrow; ++i)
		for (int j = 0; j < ncol; ++j)
			output[i][j] = input[i][j] / max_intensity * 255;
	return output;
}

bool is_interference(int i, int j, SDoublePlane &fft_real, SDoublePlane &fft_imag, double threshold) {
	double real = fft_real[i][j], imag = fft_imag[i][j];
	if (calc_log_magnitude(real, imag) > threshold)
		return true;
	return false;
}

void remove_interference_single_pixel(int i, int j, SDoublePlane &fft_real, SDoublePlane &fft_imag) {
	double threshold = 0;
	if (!is_interference(i, j, fft_real, fft_imag, threshold))  // Is not interference
		return;

	double new_real = 0, new_imag = 0;
	int pixel_count = 0;
	// Replace the pixel value by the average of the non-interfering cells in the the surrounding 8 cells
	int row_range[] = {i - 1, i - 1, i - 1,
			i,                i,
			i + 1, i + 1, i + 1};
	int col_range[] = {j - 1, j    , j + 1,
			j - 1,        j + 1,
			j - 1, j    , j + 1};

	for (int k = 0; k < 8; k++) {
		int row = row_range[k], col = col_range[k];
		if (!is_interference(row, col, fft_real, fft_imag, threshold)) {
			new_real += fft_real[row][col];
			new_imag += fft_imag[row][col];
			pixel_count++;
		}
	}
	if (pixel_count > 0) {
		fft_real[i][j] = new_real / pixel_count;
		fft_imag[i][j] = new_imag / pixel_count;
	}
}

// Write this in Part 1.2
SDoublePlane remove_interference(const SDoublePlane &input) {
	SDoublePlane fft_real, fft_imag;
	fft(input, fft_real, fft_imag);

	for (int i = 156; i <= 160; ++i)
		for (int j = 156; j <= 162; ++j)
			remove_interference_single_pixel(i, j, fft_real, fft_imag);
	for (int i = 352; i <= 356; ++i)
		for (int j = 350; j <= 356; ++j)
			remove_interference_single_pixel(i, j, fft_real, fft_imag);
	SDoublePlane output;
	ifft(fft_real, fft_imag, output);
	return output;
}

// Write this in Part 1.3 -- add watermark N to image
SDoublePlane mark_image(const SDoublePlane &input, int N) {
	int nrow = input.rows(), ncol = input.cols();
	double r = nrow / 2 * radius_fraction_CONST;  // Radius
	srandom(N);
	int v[l_CONST], topbin_x[l_CONST], topbin_y[l_CONST];
	for (int i = 0; i < l_CONST; ++i)
		v[i] = random() % 2;

	double center = nrow / 2 - 0.5;
	double theta = M_PI / (l_CONST + 1);
	for (int i = 0; i < l_CONST; ++i) {
		topbin_x[i] = (int) (center - r * cos(i * theta));
		topbin_y[i] = (int) (center + r * sin(i * theta));
	}

	SDoublePlane fft_real, fft_imag;
	fft(input, fft_real, fft_imag);

	for (int i = 0; i < l_CONST; ++i) {
		// Top half of the circle
		int row = topbin_x[i], col = topbin_y[i];
		double intensity = fft_real[row][col];
		fft_real[row][col] = intensity + alpha_CONST * v[i] * fabs(intensity);

		// Bottom half of the circle
		row = nrow - row;  // ?? computation
		col = ncol - col;
		intensity = fft_real[row][col];
		fft_real[row][col] = intensity + alpha_CONST * v[i] * fabs(intensity);
	}

	SDoublePlane output_real;
	ifft(fft_real, fft_imag, output_real);
	return output_real;
}

bool check_mark_single_pixel(const SDoublePlane &input, int row, int col, int v) {
	if (v == 0)
		return true;

	double prox_intensity = 0.125 * (
			input[row - 1][col - 1] + input[row - 1][col] + input[row - 1][col + 1] +
			input[row][col - 1] + input[row][col + 1] +
			input[row + 1][col - 1] + input[row + 1][col] + input[row + 1][col + 1]
	);
	double thr = 0.5;
	double intensity = input[row][col];
	if (fabs(intensity - prox_intensity) > thr * alpha_CONST * fabs(prox_intensity))
		return true;
	return false;
}

// Write this in Part 1.3 -- check if watermark N is in image
bool check_image(const SDoublePlane &input, int N) {
	// Re-calculate the parameters used in mark part
	int nrow = input.rows(), ncol = input.cols();
	double r = nrow / 2 * radius_fraction_CONST;  // Radius
	srandom(N);
	int v[l_CONST], topbin_x[l_CONST], topbin_y[l_CONST];
	int positive_v_count = 0;
	for (int i = 0; i < l_CONST; ++i) {
		int rand_num = random() % 2;
		v[i] = rand_num;
		positive_v_count += rand_num;
	}

	double center = nrow / 2 - 0.5;
	double theta = M_PI / (l_CONST + 1);
	for (int i = 0; i < l_CONST; ++i) {
		topbin_x[i] = (int) (center - r * cos(i * theta));
		topbin_y[i] = (int) (center + r * sin(i * theta));
	}

	SDoublePlane fft_real, fft_imag;
	fft(input, fft_real, fft_imag);

	int count = l_CONST * 2;
	for (int i = 0; i < l_CONST; ++i) {
		// Top half of the circle
		int row = topbin_x[i], col = topbin_y[i];
		if (!check_mark_single_pixel(fft_real, row, col, v[i]))
			--count;

		// Bottom half of the circle
		row = nrow - row;  // ?? computation
		col = ncol - col;
		if (!check_mark_single_pixel(fft_real, row, col, v[i]))
			--count;
	}
	printf("Detected watermark points: %d out of a total of %d\n", count - 2 * (l_CONST - positive_v_count), 2 * positive_v_count);
	return count >= 2 * (l_CONST - positive_v_count + positive_v_count * 0.5);
}


int main(int argc, char **argv)
{
	try {

		if(argc < 4)
		{
			cout << "Insufficent number of arguments; correct usage:" << endl;
			cout << "    p2 problemID inputfile outputfile" << endl;
			return -1;
		}

		string part = argv[1];
		string inputFile = argv[2];
		string outputFile = argv[3];
		cout << "In: " << inputFile <<"  Out: " << outputFile << endl;

		SDoublePlane input_image = SImageIO::read_png_file(inputFile.c_str());
		int N = atoi(argv[5]);
		/*
		 * Part 1.1
		 */
		if(part == "1.1")
		{
			// do something here!
			SDoublePlane fft_real, fft_imag;
			fft(input_image, fft_real, fft_imag);
			SDoublePlane spec_image = fft_magnitude(fft_real, fft_imag);
			if (normalize_intensity)
				spec_image = normalize_spec_img(spec_image);
			SImageIO::write_png_file(outputFile.c_str(), spec_image, spec_image, spec_image);
		}

		/*
		 * Part 1.2
		 */
		else if(part == "1.2")
		{
			// do something here!
			SDoublePlane output = remove_interference(input_image);
			SImageIO::write_png_file(outputFile.c_str(), output, output, output);
		}

		/*
		 * Part 1.3
		 */
		else if(part == "1.3")
		{
			if(argc < 6)
			{
				cout << "Need 6 parameters for watermark part:" << endl;
				cout << "    p2 1.3 inputfile outputfile operation N" << endl;
				return -1;
			}
			string op(argv[4]);
			if(op == "add")
			{
				// add watermark
				SDoublePlane output = mark_image(input_image, N);
				SImageIO::write_png_file(outputFile.c_str(), output, output, output);
				// printf("Image size: %d x %d\n", input_image.rows(), input_image.cols());
			}
			else if(op == "check")
			{
				// check watermark
				printf(check_image(input_image, N) ? ">> Watermark exists" : ">> Watermark does NOT exist");
				printf("\n===== ===== ===== ===== \n");
			}
			else
				throw string("Bad operation!");
		}
		else
			throw string("Bad part!");

	}
	catch(const string &err) {
		cerr << "Error: " << err << endl;
	}
}








