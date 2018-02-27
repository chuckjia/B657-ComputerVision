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
			kernel(i, j, 0, 0) = pixel_val;
			kernel(i, j, 0, 1) = pixel_val;
			kernel(i, j, 0, 2) = pixel_val;
		}
	return kernel;
}

void write_img_to_file(CImg<double> img, const char *filename) {
	img.normalize(0, 255).save(filename);
}

void write_level_img_to_file(CImg<double> img, const char *file_prefix, int level) {
	char filename[30];
	sprintf(filename, "z_output/%s_level_%d.jpg", file_prefix, level);
	img.normalize(0, 255).save(filename);
}

void blend() {
	int num_level = 5;
	CImg<double> M[num_level], L1[num_level], L2[num_level], G1[num_level], G2[num_level], LB[num_level];
	M[0] = CImg<double>("images/part2/mask.jpg");
	G1[0] = CImg<double>("images/part2/apple.jpg");
	G2[0]= CImg<double>("images/part2/orange.jpg");

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

	// Laplacian kernel
	double laplacian_mat[25];
	for (int i = 0; i < 25; ++i)
		laplacian_mat[i] = -gauss_mat[i];
	laplacian_mat[12] += 1;

	CImg<double> laplacian_filter = create_kernel(laplacian_mat, 5);

	for (int level = 1; level < num_level; ++level)  {
		int prev_level = level - 1;
		CImg<double> gauss_G1 = G1[level].get_convolve(gauss_filter);

		// Sub-sample
		int new_width = gauss_G1.width() / 2 - 1, new_height = gauss_G1.height() / 2 - 1;
		cimg_forXY(gauss_G1)
	}

	for (int level = 1; level < num_level; ++level) {
		// Apply filters
		prev_Mask.convolve(gauss_filter);
		prev_L1.convolve(laplacian_filter);
		prev_L2.convolve(laplacian_filter);

		// Sub-sample
		int new_width = prev_Mask.width() / 2 - 1, new_height = prev_Mask.height() / 2 - 1;
		CImg<double> Mask(new_width, new_height, 1, 3),
				L1(new_width, new_height, 1, 3),
				L2(new_width, new_height, 1, 3);

		cimg_forXY(Mask, new_x, new_y) {
			int x = new_x * 2, y = new_y * 2;
			for (int kk = 0; kk < 3; ++kk) {
				Mask(new_x, new_y, 0, kk) = prev_Mask(x, y, 0, kk);
				L1(new_x, new_y, 0, kk) = prev_L1(x, y, 0, kk);
				L2(new_x, new_y, 0, kk) = prev_L2(x, y, 0, kk);
			}
		}

		LB[level] = L1.get_mul(Mask)  - L2.get_mul(Mask - 255);
		write_level_img_to_file(Mask, "Mask", level);
		write_level_img_to_file(L1, "L1", level);
		write_level_img_to_file(L2, "L2", level);
		write_level_img_to_file(LB[level], "LB", level);
		prev_Mask = Mask;
		prev_L1 = L1;
		prev_L2 = L2;
	}

	for (int i = 0; i < 25; ++i)
		gauss_mat[i] *= 4;
	CImg<double> new_kernel = create_kernel(gauss_mat, 5);

	for (int level = num_level - 1; level > 0; --level) {
		int new_width = LB[level - 1].width(), new_height = LB[level - 1].height();
		CImg<double> LB_enlarged(new_width, new_height, 1, 3, 0);
		cimg_forXY(LB[level], x, y) {
			int new_x = 2 * x, new_y = 2 * y;
			for (int kk = 0; kk < 3; ++kk)
				LB_enlarged(new_x, new_y, 0, kk) = LB[level](x, y, 0, kk);
		}
		write_level_img_to_file(LB_enlarged, "LB_enlarged", level);
		LB[level - 1] += LB_enlarged.get_convolve(new_kernel);
		write_level_img_to_file(LB[level], "LB_test", level);
	}

	CImg<double> newMask("images/part2/mask.jpg"),
				newL1("images/part2/apple.jpg"),
				newL2("images/part2/orange.jpg");
	LB[0] = prev_L1.get_mul(newMask) - prev_L2.get_mul(newMask - 255);
	write_img_to_file(LB[0], "z_output/result.jpg");
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
