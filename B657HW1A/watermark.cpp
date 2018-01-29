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
SDoublePlane mark_image(const SDoublePlane &input, int N);

// Write this in Part 1.3 -- check if watermark N is in image
SDoublePlane check_image(const SDoublePlane &input, int N);


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

		/*
		 * Part 1.1
		 */
		if(part == "1.1")
		{
			// do something here!
			SDoublePlane fft_real, fft_imag;
			fft(input_image, fft_real, fft_imag);
			SDoublePlane spec_image = fft_magnitude(fft_real, fft_imag);
			SImageIO::write_png_file(outputFile.c_str(), spec_image, spec_image, spec_image);

			/*
			 * Testing 1: Analysis of the underlying matrix of image
			 */
			/*int nrow = spec_image.rows(), ncol = spec_image.cols();
			cout << "Image size: " << spec_image.rows() << " x " << spec_image.cols() << endl;
			FILE *image_matrix = fopen("image_matrix.txt", "wb");
			int num_high_intensity_pt = 0;
			double high_intensity_threshold = -1.5;
			for (int i = 0; i < nrow; ++i) {
				for (int j = 0; j < ncol; ++j) {
					double pixel_intensity = spec_image[i][j];
					fprintf(image_matrix, "%2.3f ", pixel_intensity);
					if (pixel_intensity > high_intensity_threshold) {
						cout << "The pixel (" << i << ", " << j << ") has high intensity " << pixel_intensity << endl;
						++num_high_intensity_pt;
					}
				}
				fprintf(image_matrix, "\n\n\n");
			}
			fclose(image_matrix);
			cout << "Total number of pixels: " << nrow * ncol
					<< " | High intensity threshold set at > " << high_intensity_threshold << endl;
			cout << "Total number of high intensity points is " << num_high_intensity_pt << endl;*/

			/*
			 * Testing 2: Find the intensity range of the interference pixels
			 */
			/*int nrow = spec_image.rows(), ncol = spec_image.cols();
			double high_intensity_threshold = 0;
			for (int i = 0; i < nrow; ++i)
				for (int j = 0; j < ncol; ++j)
					if (spec_image[i][j] < high_intensity_threshold)
						spec_image[i][j] = -10;
					else
						spec_image[i][j] = 4;
			string output_filename = "test_output.png";
			SImageIO::write_png_file(output_filename.c_str(), spec_image, spec_image, spec_image);
			printf("High intensity threshold = %1.4f\n", high_intensity_threshold);*/

			/*
			 * Testing 3: Find the intensity range of the interference pixels
			 */
			/*FILE *f = fopen("test_results/interference_area1.csv", "wb");
			int low = 150, high = 170;
			fprintf(f, "row/col,");
			for (int j = low; j < high; ++j)
				fprintf(f, "[:%d],", j);
			fprintf(f, "\n");
			for (int i = low; i < high; ++i) {
				fprintf(f, "[%d:],", i);
				for (int j = low; j < high; ++j)
					fprintf(f, "%1.4f,", spec_image[i][j]);
				fprintf(f, "\n");
			}
			fclose(f);

			f = fopen("test_results/interference_area2.csv", "wb");
			low = 340; high = 370;
			fprintf(f, "row/col,");
			for (int j = low; j < high; ++j)
				fprintf(f, "[:%d],", j);
			fprintf(f, "\n");
			for (int i = low; i < high; ++i) {
				fprintf(f, "[%d:],", i);
				for (int j = low; j < high; ++j)
					fprintf(f, "%1.4f,", spec_image[i][j]);
				fprintf(f, "\n");
			}
			fclose(f);*/
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
			}
			else if(op == "check")
			{
				// check watermark
			}
			else
				throw string("Bad operation!");

			int N = atoi(argv[5]);
		}
		else
			throw string("Bad part!");

	}
	catch(const string &err) {
		cerr << "Error: " << err << endl;
	}
}








