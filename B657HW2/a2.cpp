// B657 assignment 2 skeleton code
//
// Compile with: "make"
//
// See assignment handout for command line and project specifications.

// Author: Chuck Jia


//Link to the header file
#define cimg_use_jpeg
#include "part3.h"


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

int main(int argc, char **argv)
{
	if(argc < 6)
	{
		cout << "Insufficent number of arguments; correct usage:" << endl;
		cout << "    a2 partID image1 image2 mask num_level" << endl;
		return -1;
	}
	string part = argv[1];
	string img1_name = argv[2];
	string img2_name = argv[3];
	string mask_name = argv[4];
	int num_level = std::stoi(argv[5]);
	cout << "In: " << img1_name <<"  Out: " << img2_name << endl;

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
			func_part1("images/part1/lincoln.png");
		}
		else if(part == "part2"){
			CImg<double> left_img(img1_name.c_str());
			CImg<double> right_img(img2_name.c_str());
			CImg<double> mask_img(mask_name.c_str());
			/*CImg<double> mask_img(left_img.width(), left_img.height(), 1, 3);

			int mask_width = left_img.width();
			cimg_forXY(mask_img, x, y) {
				double left = 2 * x < mask_width ? 255 : 0;
				for (int kk = 0; kk < 3; ++kk)
					mask_img(x, y, 0, kk) = left;
			}*/

			// Blending
			blend(left_img, right_img, mask_img, num_level, true);
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
