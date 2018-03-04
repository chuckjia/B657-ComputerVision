/*
 * part3.h
 *
 *  Created on: Feb 18, 2018
 *      Author: ricciwoo
 */

#ifndef PART3_H_
#define PART3_H_
#include "part2.h"

//double threshold = 0.9;

void draw_descriptor_image(CImg<double> image,
		const vector<SiftDescriptor> descriptors, const char *filename);

int linear_interp(int x1, int x2, int y1, int y2, int x) {
	return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
}

CImg<double> merge_images(CImg<double> img1, CImg<double> img2) {
	// merge 2 images into 1
	int w1 = img1.width(), w2 = img2.width(), w = w1 + w2;
	int h = img1.height(), c = img1.spectrum();
	CImg<double> image(w, h, 1, 3);

	for (int kk = 0; kk < c; ++kk)
		for (int j = 0; j < h; ++j) {
			for (int i = 0; i < w1; ++i)
				image(i, j, 0, kk) = img1(i, j, 0, kk);
			for (int i = 0; i < w2; ++i)
				image(i + w1, j, 0, kk) = img2(i, j, 0, kk);
		}
	return image;
}

void draw_match_lines(CImg<double> img1, CImg<double> img2,
		vector<pair<SiftDescriptor, SiftDescriptor> > &descr, string filename) {

	// merge 2 images into 1
	CImg<double> image = merge_images(img1, img2);

	// draw a line for each match pair
	int n = descr.size(), width = img1.width();
	double color_point[] = {255.0, 255.0, 0};
	for (int d_no = 0; d_no < n; ++d_no) {
		int c1 = descr[d_no].first.col, c2 = width + descr[d_no].second.col;
		int r1 = descr[d_no].first.row, r2 = descr[d_no].second.row;
		for (int x = c1; x <= c2; ++x) {
			int y = linear_interp(c1, c2, r1, r2, x);
			for (int kk = 0; kk < 3; ++kk)
				image(x, y, 0, kk) = color_point[kk];
		}
	}
	image.save(filename.c_str());
}

void calc_match_pairs(vector<SiftDescriptor> &descr_src, vector<SiftDescriptor> &descr_dst,
		vector<pair<SiftDescriptor, SiftDescriptor> > &descr_match, double threshold) {
	// Calculate matches (closest / second-closest < threshold)
	int dim_src = descr_src.size(), dim_dst = descr_dst.size();
	for (int i = 0; i < dim_src; ++i) {
		int match_index = 0;
		double closest = 1e10, second = closest;

		SiftDescriptor descr_i = descr_src[i];
		for (int j = 0; j < dim_dst; ++j) {
			double dist = 0.0;
			SiftDescriptor descr_j = descr_dst[j];
			for (int k = 0; k < 128; ++k)
				dist += pow(descr_i.descriptor[k] - descr_j.descriptor[k], 2);
			dist = sqrt(dist);

			if (dist < closest) {
				second = closest;
				closest = dist;
				match_index = j;
			} else if (dist < second) {
				second = dist;
			}
		}
		if (closest / second < threshold)
			descr_match.push_back(pair<SiftDescriptor, SiftDescriptor>(descr_i,
					descr_dst[match_index]));
	}
	printf("descr_match.size() = %d\n", descr_match.size());
}

void ransac(vector<pair<SiftDescriptor, SiftDescriptor> > &descr, double diff) {
	// simple implementation, use each consecutive 4 pairs in descriptors
	int n = descr.size(), m = n - 3;
	double homography[9];
	bool inliners[n];
	memset(inliners, 0, sizeof(bool) * n);
	int i_sel = 0;
	int max_count = 0;
	for (int i = 0; i < m; ++i) {
		// calculate homography
		double pointsA[4][2], pointsB[4][2];
		for (int j = 0; j < 4; ++j) {
			pointsA[j][0] = descr[i + j].first.col;
			pointsA[j][1] = descr[i + j].first.row;
			pointsB[j][0] = descr[i + j].second.col;
			pointsB[j][1] = descr[i + j].second.row;
		}
		solving_homography(pointsA, pointsB, (double*)homography);
		// calculate votes
		int count = 0;
		for (int j = 0; j < n; ++j) {
			double tempt[3], result[3];
			tempt[0] = descr[j].first.col;
			tempt[1] = descr[j].first.row;
			tempt[2] = 1.0;
			for (int k = 0; k < 3; ++k) {
				result[k] = 0.0;
				for (int l = 0; l < 3; ++l) {
					result[k] += homography[k*3 + l] * tempt[l];
				}
			}
			result[0] /= result[2];
			result[1] /= result[2];
			if (sqrt(pow(result[0] - descr[j].second.col, 2.0)
					+ pow(result[1] - descr[j].second.row, 2.0)) < diff) {
				count ++;
			}
		}
//		int minCount = n / 2;
//		if (count > minCount) {
//			for (int j = 0; j < 4; ++j) {
//				inliners[i + j] = true;
//			}
//		}
		if (count > max_count) {
			max_count = count;
			i_sel = i;
		}
	}
	// calculate inliners
	vector<pair<SiftDescriptor, SiftDescriptor> > descr_ransac;
//	for (int i = 0; i < n; ++i) {
//		if (inliners[i]) {
//			descr_ransac.push_back(descr[i]);
//		}
//	}
	for (int i = 0; i < 4; ++i) {
		descr_ransac.push_back(descr[i_sel + i]);
	}
	// update descriptor
	int s = descr_ransac.size();
	descr.clear();
	for (int i = 0; i < s; ++i) {
		descr.push_back(descr_ransac[i]);
	}
}

void func_part3(string inFile_src, string inFile_dst, double threshold, double diff){

	// calculate feature points of image_src, using SIFT
	CImg<double> image_src(inFile_src.c_str());
	CImg<double> gray_src = image_src.get_RGBtoHSI().get_channel(2);
	vector<SiftDescriptor> descr_src = Sift::compute_sift(gray_src);
	//draw_descriptor_image(image_src, descr_src, "descr_src.png");

	// calculate feature points of image_dst, using SIFT
	CImg<double> image_dst(inFile_dst.c_str());
	CImg<double> gray_dst = image_dst.get_RGBtoHSI().get_channel(2);
	vector<SiftDescriptor> descr_dst = Sift::compute_sift(gray_dst);
	//draw_descriptor_image(image_dst, descr_dst, "descr_dst.png");

	// calculate match pairs
	vector<pair<SiftDescriptor, SiftDescriptor> > descr_match;
	calc_match_pairs(descr_src, descr_dst, descr_match, threshold);

	// output matches pairs from SIFT
	string filename = "sift_matches.jpg";
	draw_match_lines(image_src, image_dst, descr_match, filename);

	// RANSAC
	ransac(descr_match, diff);

	// output matches pairs from RANSAC
	filename = "ransc_matches.jpg";
	draw_match_lines(image_src, image_dst, descr_match, filename);
}

#endif /* PART3_H_ */
