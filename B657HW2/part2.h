/*
 * part2.h
 *
 *  Created on: Feb 20, 2018
 *      Author: Chuck Jia
 */

#ifndef PART2_H_
#define PART2_H_

#include "part1.h"

// Create filter kerner with an ndim x ndim matrix, as an 1D array
CImg<double> create_kernel(double mat[], int ndim) {
	CImg<double> kernel(ndim, ndim, 1, 3);
	for (int i = 0; i < ndim; ++i)
		for (int j = 0; j < ndim; ++j) {
			double pixel_val = mat[i * ndim + j];
			for (int kk = 0; kk < 3; ++kk)
				kernel(i, j, 0, kk) = pixel_val;
		}
	return kernel;
}

// Write CImg object to image file
void write_img_to_file(CImg<double> &img, const char *filename) {
	img.get_normalize(0, 255).save(filename);
}

// Write CImg object to image file, without normalization
void write_img_to_file_unnormalized(CImg<double> &img, const char *filename) {
	img.save(filename);
}

// Write one level of a pyramid to image file
void write_level_img_to_file(CImg<double> &img, const char *file_prefix, int level) {
	char filename[30];
	sprintf(filename, "z_output/%s_level_%d.jpg", file_prefix, level);
	img.get_normalize(0, 255).save(filename);
}

// Write the matrix for an CImg object to csv file
void write_img_mat_to_file(CImg<double> &img, const char *file_prefix) {
	char r_name[30], g_name[30], b_name[30];
	sprintf(r_name, "z_output/%s_R.csv", file_prefix);
	sprintf(g_name, "z_output/%s_G.csv", file_prefix);
	sprintf(b_name, "z_output/%s_B.csv", file_prefix);
	FILE *rf = fopen(r_name, "wb");
	FILE *gf = fopen(g_name, "wb");
	FILE *bf = fopen(b_name, "wb");
	int width = img.width(), height = img.height();
	for (int y = 0; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
			fprintf(rf, "%3.4f,", img(x, y, 0, 0));
			fprintf(gf, "%3.4f,", img(x, y, 0, 1));
			fprintf(bf, "%3.4f,", img(x, y, 0, 2));
		}
		fprintf(rf, "\n");
		fprintf(gf, "\n");
		fprintf(bf, "\n");
	}
	fclose(rf); fclose(gf); fclose(bf);
}

// Blend two images using a mask image
// Input: left_img, right_img are the two images to blend
//        mask_img is the mask image. Its pixel values are in {0, 255} with 0 representing one image and 255 the other
//        print_result is a bool, indicating if the pyramid images and the result images are to be printed as files
// Output: The blended image as a CImg object
CImg<double> blend(CImg<double> &left_img, CImg<double> &right_img, CImg<double> mask_img, bool print_result) {
	/*
	 * Create filters need in blending
	 */

	// Matrix for the 5 x 5 Gaussian kernel
	double gauss_mat[] = { 1, 4, 6, 4, 1,
			4, 16, 24, 16, 4,
			6, 24, 36, 24, 6,
			4, 16, 24, 16, 4,
			1,  4,  6,  4, 1 };
	for (int i = 0; i < 25; ++i) gauss_mat[i] *= 1 / 256.;
	CImg<double> gauss_filter = create_kernel(gauss_mat, 5);

	// Four times Gaussian kernel. Used to approximate for upscaled images
	for (int i = 0; i < 25; ++i) gauss_mat[i] *= 4;
	CImg<double> gauss_upscale_kern = create_kernel(gauss_mat, 5);

	/*
	 * Pyramid blending
	 */

	int num_level = 5;  // Number of levels used in pyramids
	CImg<double> M[num_level], G1[num_level], G2[num_level];  // M, G1, G2 are Gaussian pyramids for the mask, image 1 and image 2

	/*
	 * Create Gaussian pyramids
	 */

	G1[0] = left_img; G2[0]= right_img; M[0] = mask_img;

	for (int level = 1; level < num_level; ++level)  {
		int prev_level = level - 1;
		CImg<double> gauss_M = M[prev_level].get_convolve(gauss_filter, false, false);
		CImg<double> gauss_G1 = G1[prev_level].get_convolve(gauss_filter, false, false);
		CImg<double> gauss_G2 = G2[prev_level].get_convolve(gauss_filter, false, false);

		// Sub-sample
		int orig_width = gauss_M.width(), orig_height = gauss_M.height();
		int new_width = orig_width % 2 ? orig_width / 2 + 1 : orig_width / 2,
				new_height = orig_height % 2 ? orig_height / 2 + 1 : orig_height / 2;
		M[level] = CImg<double>(new_width, new_height, 1, 3);
		G1[level] = CImg<double>(new_width, new_height, 1, 3);
		G2[level] = CImg<double>(new_width, new_height, 1, 3);

		cimg_forXY(M[level], new_x, new_y) {
			int x = new_x * 2, y = new_y * 2;
			for (int kk = 0; kk < 3; ++kk) {
				M[level](new_x, new_y, 0, kk) = gauss_M(x, y, 0, kk);
				G1[level](new_x, new_y, 0, kk) = gauss_G1(x, y, 0, kk);
				G2[level](new_x, new_y, 0, kk) = gauss_G2(x, y, 0, kk);
			}
		}
	}

	// Write the Gaussian pyramids to file
	if (print_result)
		for (int level = 0; level < num_level; ++level) {
			write_level_img_to_file(M[level], "M", level);
			write_level_img_to_file(G1[level], "G1", level);
			write_level_img_to_file(G2[level], "G2", level);
		}

	/*
	 * Create Laplacian pyramids
	 */

	CImg<double> L1[num_level], L2[num_level];  // L1 and L2 are the Laplacian pyramids for image 1 and image 2
	int last_level = num_level - 1;
	for (int level = 0; level < last_level; ++level) {
		L1[level] = G1[level] - G1[level].get_convolve(gauss_filter, false, false);
		L2[level] = G2[level] - G2[level].get_convolve(gauss_filter, false, false);
	}

	L1[last_level] = G1[last_level];
	L2[last_level] = G2[last_level];

	// Write the Laplacian pyramids to file
	if (print_result)
		for (int level = 0; level < num_level; ++level) {
			write_level_img_to_file(L1[level], "L1", level);
			write_level_img_to_file(L2[level], "L2", level);
		}

	/*
	 *	Creating "blended" Laplacian pyramid
	 */

	CImg<double> LB[num_level];  // LB is the "Laplacian" pyramid for the blended image
	for (int level = 0; level < num_level; ++level) {
		LB[level] = CImg<double>(M[level].width(), M[level].height(), 1, 3);
		cimg_forXY(M[level], x, y) {
			for (int kk = 0; kk < 3; ++kk) {
				double M_val = M[level](x, y, 0, kk) / 255.;
				LB[level](x, y, 0, kk) = L1[level](x, y, 0, kk) * M_val + L2[level](x, y, 0, kk) * (1. - M_val);
			}
		}
	}

	// Write the pyramid to file
	if (print_result)
		for (int level = 0; level < num_level; ++level)
			write_level_img_to_file(LB[level], "LB", level);

	/*
	 *	Combine all images
	 */
	for (int level = num_level - 1; level > 0; --level) {
		printf("Progress: level no. %d\n", level);
		int prev_level = level - 1, new_width = LB[prev_level].width(), new_height = LB[prev_level].height();

		// Upscaling
		CImg<double> new_img(new_width, new_height, 1, 3, 0);
		cimg_forXY(LB[level], x, y) {
			int new_x = 2 * x, new_y = 2 * y;
			for (int kk = 0; kk < 3; ++kk)
				new_img(new_x, new_y, 0, kk) = LB[level](x, y, 0, kk);
		}

		new_img.convolve(gauss_upscale_kern, false, false);
		LB[prev_level] += new_img;

		// write_level_img_to_file(LB[prev_level], "Final_LB", prev_level);
	}

	if (print_result)
		write_img_to_file_unnormalized(LB[0], "z_output/Final_Result.jpg");
	return LB[0];
}

#endif /* PART2_H_ */
