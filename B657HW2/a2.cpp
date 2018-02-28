// B657 assignment 2 skeleton code
//
// Compile with: "make"
//
// See assignment handout for command line and project specifications.

// Author: Chuck Jia


//Link to the header file
#define cimg_use_jpeg
#include "CImg.h"
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Sift.h>


//Use the cimg namespace to access functions easily
using namespace cimg_library;
using namespace std;

void draw_descriptor_image(CImg<double> image, const vector<SiftDescriptor> descriptors, const char *filename)
{
	for(unsigned int i=0; i < descriptors.size(); i++)
	{
		int tx1 = 0, ty1 = 0, tx2 = 0, ty2 = 0;
		double color_point[] = {255.0, 255.0, 0};
		for(int x=-2; x<3; x++)
			for(int y=-2; y<3; y++)
				if(x==0 || y==0)
					for(int c=0; c<3; c++){
						//Find if coordinates are in workspace to draw crosshair
						tx1 = (descriptors[i].col + y - 1);
						ty1 = (descriptors[i].row + x - 1);
						if (tx1 >= 0 && tx1 < image.width() && ty1 >= 0 && ty1 < image.height())
							image( tx1, ty1, 0, c) = color_point[c];
					}
	}
	image.get_normalize(0,255).save(filename);
}

void testmsg() {
	printf("\n===== ===== ===== ===== ===== \n>> This is a test message.\n===== ===== ===== ===== ===== \n\n");
}

// Author: Chuck Jia
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

void write_img_to_file(CImg<double> &img, const char *filename) {
	img.normalize(0, 255).save(filename);
}

void write_level_img_to_file(CImg<double> &img, const char *file_prefix, int level) {
	char filename[30];
	sprintf(filename, "z_output/%s_level_%d.jpg", file_prefix, level);
	img.normalize(0, 255).save(filename);
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
	for (int x = 0; x < width; ++x) {
		for (int y = 0; y < height; ++y) {
			fprintf(rf, "%3.4f ", img(x, y, 0, 0));
			fprintf(gf, "%3.4f ", img(x, y, 0, 1));
			fprintf(bf, "%3.4f ", img(x, y, 0, 2));
		}
		fprintf(rf, "\n");
		fprintf(gf, "\n");
		fprintf(bf, "\n");
	}
	fclose(rf); fclose(gf); fclose(bf);
}

void blend() {

	/*cimg_forXY(prev_Mask, x, y) {
		// printf("(%d, %d) = (%f, %f, %f)\n",
		//x, y, prev_Mask(x, y, 0, 0), prev_Mask(x, y, 0, 1), prev_Mask(x, y, 0, 2));
		for (int kk = 0; kk < 3; ++kk)
			prev_Mask(x, y, 0, kk) = 2 * x < prev_Mask.width() ? 0 : 1;
	}*/

	// Gaussian kernel
	double gauss_mat[] = {1, 4, 6, 4, 1,
			4, 16, 24, 16, 4,
			6, 24, 36, 24, 6,
			4, 16, 24, 16, 4,
			1,  4,  6,  4, 1};
	for (int i = 0; i < 25; ++i)
		gauss_mat[i] *= 1 / 256.;
	CImg<double> gauss_filter = create_kernel(gauss_mat, 5);

	// write_img_mat_to_file(gauss_filter, "gauss_filter");

	// Laplacian kernel
	double laplacian_mat[25];
	for (int i = 0; i < 25; ++i)
		laplacian_mat[i] = -gauss_mat[i];
	laplacian_mat[12] += 1;

	CImg<double> laplacian_filter = create_kernel(laplacian_mat, 5);


	int num_level = 6;
	CImg<double> M[num_level], L1[num_level], L2[num_level], G1[num_level], G2[num_level], LB[num_level];


	/*
	 * Creating Gaussian pyramids
	 */

	// M[0] = CImg<double>("images/part2/mask.jpg");
	G1[0] = CImg<double>("images/part2/apple.jpg");
	G2[0]= CImg<double>("images/part2/orange.jpg");
	M[0] = CImg<double>(G1[0].width(), G1[0].height(), 1, 3);

	int orig_width = G1[0].width();

	cimg_forXY(M[0], x, y) {
		int right = (2 * x < orig_width) * 255;
		for (int kk = 0; kk < 3; ++kk)
			M[0](x, y, 0, kk) = right;
	}
	write_level_img_to_file(M[0], "M", 0);

	for (int level = 1; level < num_level; ++level)  {
		int prev_level = level - 1;
		CImg<double> gauss_M = M[prev_level].get_convolve(gauss_filter);
		CImg<double> gauss_G1 = G1[prev_level].get_convolve(gauss_filter);
		CImg<double> gauss_G2 = G2[prev_level].get_convolve(gauss_filter);

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
		write_level_img_to_file(M[level], "M", level);
		write_level_img_to_file(G1[level], "G1", level);
		write_level_img_to_file(G2[level], "G2", level);
	}

	/*
	 * Creating Laplacian pyramids
	 */

	int last_level = num_level - 1;
	for (int level = 0; level < last_level; ++level)  {
		int new_width = G1[level].width(), new_height = G1[level].height();
		L1[level] = G1[level] - G1[level].get_convolve(gauss_filter);
		L2[level] = G2[level] - G2[level].get_convolve(gauss_filter);

		// int next_level = level + 1;
		/*cimg_forXY(G1[level], x, y) {
			bool not_from_upscale = (x % 2) || (y % 2);
			for (int kk = 0; kk < 3; ++kk) {
				L1[level](x, y, 0, kk) = G1[level](x, y, 0, kk) - (not_from_upscale ? 0 : G1[next_level](x / 2, y / 2, 0, kk));
				L2[level](x, y, 0, kk) = G2[level](x, y, 0, kk) - (not_from_upscale ? 0 : G2[next_level](x / 2, y / 2, 0, kk));
			}
		}*/
		write_level_img_to_file(L1[level], "L1", level);
		write_level_img_to_file(L2[level], "L2", level);
	}

	L1[last_level] = G1[last_level];
	L2[last_level] = G2[last_level];

	/*
	 *	Creating "blended" Laplacian pyramid
	 */

	for (int level = 0; level < num_level; ++level) {
		LB[level] = CImg<double>(M[level].width(), M[level].height(), 1, 3);
		// CImg<double> left = L1[level].get_mul(M[level]), right = L2[level].get_mul(1 - M[level]);
		cimg_forXY(M[level], x, y) {
			for (int kk = 0; kk < 3; ++kk) {
				double M_val = M[level](x, y, 0, kk);
				LB[level](x, y, 0, kk) = L1[level](x, y, 0, kk) * M_val / 255. + L2[level](x, y, 0, kk) * (255 - M_val) / 255.;
			}
		}
		// LB[level] = L1[level].get_mul(M[level]) + L2[level].get_mul(1 - M[level]);
		//write_level_img_to_file(left, "left", level);
		//write_level_img_to_file(right, "right", level);
		write_level_img_to_file(LB[level], "LB", level);
	}

	/*
	 *	Combining all images
	 */
	for (int i = 0; i < 25; ++i)
		gauss_mat[i] *= 4;
	CImg<double> new_kernel = create_kernel(gauss_mat, 5);

	for (int level = num_level - 1; level > 0; --level) {
		printf("Level no. %d", level);
		int prev_level = level - 1;
		int new_width = LB[prev_level].width(), new_height = LB[prev_level].height();

		// Upscaling
		CImg<double> new_img(new_width, new_height, 1, 3, 0);
		cimg_forXY(LB[level], x, y) {
			int new_x = 2 * x, new_y = 2 * y;
			for (int kk = 0; kk < 3; ++kk)
				new_img(new_x, new_y, 0, kk) = LB[level](x, y, 0, kk);
		}
		new_img.convolve(gauss_filter);

		// LB[prev_level] += LB[level].get_resize(LB[prev_level], 3);
		LB[prev_level] += new_img;
		write_level_img_to_file(LB[prev_level], "Final_LB", prev_level);
	}


}

int main()//int argc, char **argv)
{
	try {

		/*
      	  TEST CODE - STARTS
		 */
		string part = "part2";
		/*CImg<double> input_image("images/part2/apple.jpg");
		CImg<double> input_gray = input_image.get_RGBtoHSI().get_channel(2);
		vector<SiftDescriptor> input_descriptors = Sift::compute_sift(input_gray);
		draw_descriptor_image(input_image, input_descriptors, "input_image.jpg");*/
		/*
      	  TEST CODE - ENDS
		 */

		if(part == "part1"){
			// Billboard
		}
		else if(part == "part2"){
			// Blending
			blend();
		}
		else if(part == "part3"){
			// RANSAC
		}
		else if(part == "part4"){
			// Panorama
		}


		// feel free to add more conditions for other parts (e.g. more specific)
		//  parts, for debugging, etc.
	}
	catch(const string &err) {
		cerr << "Error: " << err << endl;
	}
}
