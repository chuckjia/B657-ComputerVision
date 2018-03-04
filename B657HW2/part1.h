/*
 * part1.h
 *
 *  Created on: Feb 20, 2018
 *      Author: ricciwoo
 */

#ifndef PART1_H_
#define PART1_H_
#include "CImg.h"
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Sift.h>

// Matrix multiplication
void matrix_multiplication(double A[3][3], double B[3], double result[3]) {
	for (int i = 0; i < 3; i++) {
		result[i] = 0;
		for (int j = 0; j < 3; j++)
			result[i] += A[i][j] * B[j];
	}
}

double TansformMatrixInvert[3][3] = { {1.1246685805, -0.3146766040, 222.940925},
		{0.1088390506, 0.6850586647, -19.924695},
		{0.0002645872, -0.0005970689, 1.082785} };

// Transformation matrix calculation
void solving_homography(double pointsA[4][2], double pointsB[4][2], double* homography) {
	// Cramer's rule
	// ax + by + cz = j          | a b c |
	// dx + ey + fz = k    det = | d e f |
	// gx + hy + iz = l          | g h i |
	//     | j b c |            | a j c |
	// x = | k e f | / det, y = | d k f | / det, z = ...
	//     | l h i |            | g l i |
	CImg<double> coeff(8, 8, 1, 1);
	double right[8], tempt[8], result[8];
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 2; ++j) {
			coeff(j,   2*i,   0, 0) = pointsA[i][j];
			coeff(j,   2*i+1, 0, 0) = 0.0;
			coeff(3+j, 2*i,   0, 0) = 0.0;
			coeff(3+j, 2*i+1, 0, 0) = pointsA[i][j];
			coeff(6+j, 2*i,   0, 0) = - pointsA[i][j]*pointsB[i][0];
			coeff(6+j, 2*i+1, 0, 0) = - pointsA[i][j]*pointsB[i][1];
		}
		coeff(2, 2*i,   0, 0) = 1.0;
		coeff(2, 2*i+1, 0, 0) = 0.0;
		coeff(5, 2*i,   0, 0) = 0.0;
		coeff(5, 2*i+1, 0, 0) = 1.0;
		right[2*i]   = pointsB[i][0];
		right[2*i+1] = pointsB[i][1];
	}
	// write matrix form of each pair of point (x0,y0) and (x0',y0')
	// note that must write (x0'w0, y0'w0, w0) for homogeneous coord
	// a b c     x0     x0'w0       ax0 + by0 + c = x0'w0
	// d e f  *  y0  =  y0'w0  ==>  dx0 + ey0 + f = y0'w0  ==>
	// g h 1     1       w0         gx0 + hy0 + 1 = w0
	// -----------------------------------------
	// ax0 + by0 + c - gx0x0' - hy0x0' = x0'
	// dx0 + ey0 + f - gx0y0' + hy0y0' = y0' ==>
	// -----------------------------------------
	// x0a + y0b + 1c +  0d +  0e + 0f - x0x0'g - y0x0'h = x0'
	//  0a +  0b + 0c + x0d + y0e + 1f - x0y0'g - y0y0'h = y0'
	// same for (xi,yi) and (xi',yi'), i = 1, 2, 3
	// -------------------------------------------
	// Solve for 8 dimensional linear system
	double deter = coeff.det();
	for (int i = 0; i < 8; ++i) {
		for (int j = 0; j < 8; ++j) {
			tempt[j] = coeff(i, j, 0, 0);
			coeff(i, j, 0, 0) = right[j];
		}
		result[i] = coeff.det() / deter;
		for (int j = 0; j < 8; ++j) {
			coeff(i, j, 0, 0) = tempt[j];
		}
	}
	// set output array
	for (int i = 0; i < 8; ++i) {
		homography[i] = result[i];
	}
	homography[8] = 1.0;
}

//Inverse Warpping
CImg<double> inverseWarpping(CImg<double> orig_img, CImg<double> new_img, double orig_coor[4][2], double new_coor[4][2]) {
	double homography[9];
	solving_homography(new_coor, orig_coor, (double*)homography);
	double Trans[3][3];
	int index = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Trans[i][j] = homography[index++];
		}
	}
	cimg_forXYC(new_img,x,y,c) {
		double orig_coor[] = {x,y,1};
		double new_coor[3];
		matrix_multiplication(Trans, orig_coor, new_coor);
		double n_x = new_coor[0]/new_coor[2];
		double n_y = new_coor[1]/new_coor[2];
		if (n_x >= 0 && n_x < orig_img.width() && n_y >= 0 && n_y < orig_img.height()) {
			new_img(x, y, 0, c) = orig_img(n_x, n_y, 0, c);
		}
	}
	return new_img;
}

void func_part1(string poster_img) {
	//1.1
	//	//This part input command example: "./a2 part1 lincoln.png"

	//1.2
	CImg<double> book2("images/part1/book2.jpg");
	CImg<double> book2_warped(book2.width(), book2.height(), 1, 3);

	//Caluculating transform matrix
	double pointsA[4][2] = {{318, 256}, {534, 372}, {316, 670}, {73, 473}};
	double pointsB[4][2] = {{141, 131}, {480, 159}, {493, 630}, {64, 601}};
	book2_warped = inverseWarpping(book2, book2_warped, pointsB, pointsA);
	book2_warped.save("book2_warpped.jpg");

	//1.3
	//input image information:
	CImg<double> part1_3_input(poster_img.c_str());
	int input_width = part1_3_input.width();
	int input_height = part1_3_input.height();
	double part1_3_input_coor[4][2] = {{0, 0}, {input_width, 0},
			{0, input_height}, {input_width, input_height}};

	//billboard1 result
	CImg<double> billboard1("images/part1/billboard1.jpg");
	double part1_3_billboard1_coor[4][2] = {{102, 61}, {532, 61}, {102, 203}, {532, 203}};
	billboard1 = inverseWarpping(part1_3_input, billboard1, part1_3_input_coor,
			part1_3_billboard1_coor);
	billboard1.save("newbillboard1.jpg");

	//billboard3 result
	CImg<double> billboard3("images/part1/billboard3.jpg");
	double part1_3_billboard3_coor[4][2] = {{616, 287}, {1259, 261}, {610, 608}, {1261, 604}};
	billboard3 = inverseWarpping(part1_3_input, billboard3, part1_3_input_coor,
			part1_3_billboard3_coor);
	billboard3.save("newbillboard3.jpg");
}

#endif /* PART1_H_ */
