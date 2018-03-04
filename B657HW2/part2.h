/*
 * part2.h
 *
 *  Created on: Feb 20, 2018
 *      Author: chuckjia
 */

#ifndef PART2_H_
#define PART2_H_

#include "part1.h"

void testmsg() {
	printf("\n===== ===== ===== ===== ===== \n>> This is a test message.\n===== ===== ===== ===== ===== \n\n");
}

// Author: Chuck Jia
CImg<double> create_kernel(double mat[5][5]) {
	int ndim = 5;
	CImg<double> kernel(ndim, ndim, 1, 3);
	for (int i = 0; i < ndim; ++i)
		for (int j = 0; j < ndim; ++j) {
			double pixel_val = mat[i][j];
			for (int kk = 0; kk < 3; ++kk)
				kernel(i, j, 0, kk) = pixel_val;
		}
	return kernel;
}

void write_img_to_file(CImg<double> &img, const char *filename) {
	img.get_normalize(0, 255).save(filename);
}

void write_img_to_file_unnormalized(CImg<double> &img, const char *filename) {
	img.save(filename);
}

void write_level_img_to_file(CImg<double> &img, const char *file_prefix, int level) {
	char filename[30];
	sprintf(filename, "z_output/%s_level_%d.jpg", file_prefix, level);
	img.get_normalize(0, 255).save(filename);
}

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

CImg<double> calc_corr(CImg<double> &image, double kernel[5][5]) {
	int width = image.width(), height = image.height();
	CImg<double> ans(width, height, 1, 3);
	cimg_forXY(image, x, y) {
		for (int kk = 0; kk < 3; ++kk) {
			double val = 0;
			for (int i = -2; i <= 2; ++i) {
				int ii = i + 2, xx = x + i;
				xx = xx < 0 ? -xx : xx;
				xx = xx >= width ? 2 * width - 2 - xx : xx;
				for (int j = -2; j <= 2; ++j) {
					int jj = j + 2, yy = y + j;
					yy = yy < 0 ? -yy : yy;
					yy = yy >= height ? 2 * height - 2 - yy : yy;
					val += image(xx, yy, 0, kk) * kernel[ii][jj];
				}
			}
			ans(x, y, 0, kk) = val;
		}
	}
	return ans;
}

void blend(CImg<double> &left_img, CImg<double> &right_img, CImg<double> mask_img) {
	/*
	 * Create filters need in blending
	 */

	// Matrix for the 5 x 5 Gaussian kernel
	double gauss_mat[5][5] = { {1, 4, 6, 4, 1},
			{4, 16, 24, 16, 4},
			{6, 24, 36, 24, 6},
			{4, 16, 24, 16, 4},
			{1,  4,  6,  4, 1} };
	for (int i = 0; i < 5; ++i)
		for (int j = 0; j < 5; ++j)
			gauss_mat[i][j] *= 1 / 256.;
	CImg<double> gauss_filter = create_kernel(gauss_mat);

	// Four times Gaussian kernel. Used to approximate for upscaled images
	double gauss_times_4_mat[5][5];
	for (int i = 0; i < 5; ++i)
		for (int j = 0; j < 5; ++j)
			gauss_times_4_mat[i][j] = gauss_mat[i][j] * 4;
	CImg<double> new_kernel = create_kernel(gauss_times_4_mat);

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
		int new_width = gauss_M.width() / 2, new_height = gauss_M.height() / 2;
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
	for (int level = 0; level < num_level; ++level) {
		write_level_img_to_file(M[level], "M", level);
		write_level_img_to_file(G1[level], "G1", level);
		write_level_img_to_file(G2[level], "G2", level);
	}

	/*
	 * Creating Laplacian pyramids
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
	for (int level = 0; level < num_level; ++level)
		write_level_img_to_file(LB[level], "LB", level);

	/*
	 *	Combine all images
	 */
	for (int level = num_level - 1; level > 0; --level) {
		printf("Level no. %d\n", level);
		int prev_level = level - 1;
		// int new_width = LB[prev_level].width(), new_height = LB[prev_level].height();
		int new_width = 2 * LB[level].width(), new_height = 2 * LB[level].height();

		// Upscaling
		CImg<double> new_img(new_width, new_height, 1, 3, 0);
		cimg_forXY(LB[level], x, y) {
			int new_x = 2 * x, new_y = 2 * y;
			for (int kk = 0; kk < 3; ++kk)
				new_img(new_x, new_y, 0, kk) = LB[level](x, y, 0, kk);
		}

		/*if (new_width % 2) {
			int last_col = new_width - 1;
			for (int y = 0; y < new_height; ++y)
				if (y % 2 == 1)
					for (int kk = 0; kk < 3; ++kk)
						new_img(last_col, y, 0, kk) = G2[prev_level](last_col, y, 0, kk);
		}*/

		// new_img.convolve(new_kernel, false, false);
		new_img = calc_corr(new_img, gauss_times_4_mat);
		// LB[prev_level] += new_img;
		cimg_forXY(new_img, x, y) {
			for (int kk = 0; kk < 3; ++kk) {
				LB[prev_level](x, y, 0, kk) += new_img(x, y, 0, kk);
			}
		}

		write_level_img_to_file(LB[prev_level], "Final_LB", prev_level);
	}

	write_img_to_file_unnormalized(LB[0], "z_output/Final_Result.jpg");
}

#endif /* PART2_H_ */
