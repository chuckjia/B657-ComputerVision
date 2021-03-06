#! /usr/bin/env

#############################################################
# B657 Homework Assignment #1 IC Detection Evaluation Script
# Usage: -> python a1-eval.py gt.txt detected.txt
#############################################################

import sys;
import numpy as np;
from matplotlib import pyplot as pl;

if len(sys.argv) != 3:
	print "Usage: -> python a1-eval.py gt.txt detected.txt";
	sys.exit();

# read ground truth file
gt_file = open(sys.argv[1], 'r');
gt_lines = gt_file.readlines();
gt_file.close();

# read detection file
dt_file = open(sys.argv[2], 'r');
dt_lines = dt_file.readlines();
dt_file.close();

# extract groundtruths & detections
dtype = [('row', int), ('col', int), ('height', int), ('width', int), ('confidence', float)];
gt_array = np.zeros((len(gt_lines),),dtype=dtype);
dt_array = np.zeros((len(dt_lines),),dtype=dtype);

cnt = 0;
for line in gt_lines:
	items = line.rsplit(" ");
	gt_array[cnt] = (int(items[0]), int(items[1]), int(items[2]), int(items[3]), float(items[4]));
	cnt = cnt + 1;

cnt = 0;
for line in dt_lines:
        items = line.rsplit(" ");
        dt_array[cnt] = (int(items[0]), int(items[1]), int(items[2]), int(items[3]), float(items[4]));
        cnt = cnt + 1;

# sort detection array
sorted_dt_array = np.sort(dt_array, None, 'quicksort', 'confidence')[::-1];

# define symbol set
symbol_set = ['ic'];

# use top K detections to evaluate
ov_thresh = 0.5;

# mean average precision (mAP)
maps = (0, 0, 0);
mAP = 0.0;

# total TP and #correct pitches (only calculated when recall is maximum)
total_tp = 0;
num_gt_symbol_type = 0;

for i in range(0,len(symbol_set)):
	symbol = symbol_set[i];
	sorted_dt_array_symbol = sorted_dt_array;
	gt_array_symbol = gt_array;
        
        num_gt_symbol_type = num_gt_symbol_type + 1;
	pr_array = np.zeros((len(sorted_dt_array_symbol),),dtype=float);
	rc_array = np.zeros((len(sorted_dt_array_symbol),),dtype=float);
	# assign detection to ground truth
	dt_assign = np.zeros((len(sorted_dt_array_symbol),),dtype=('int'));
	# initialize with -1
	for k in range(0,len(dt_assign)):
		dt_assign[k] = -1;
	gt_cnt = 0;
	for gt in gt_array_symbol:
		gt_cnt = gt_cnt + 1;
		dt_cnt = 0;
		detected = False;
		for dt in sorted_dt_array_symbol:
			dt_cnt = dt_cnt + 1;
			# skip detections already assigned
			if dt_assign[dt_cnt - 1] != -1:
				continue;
			# calculate intersection-over-union ratio (IOU)
			xx1i = max(gt['col'], dt['col']);
			yy1i = max(gt['row'], dt['row']);
			xx1u = min(gt['col'], dt['col']);
			yy1u = min(gt['row'], dt['row']);
			xx2i = min(gt['col'] + gt['width'], dt['col'] + dt['width']);
			yy2i = min(gt['row'] + gt['height'], dt['row'] + gt['height']);
			xx2u = max(gt['col'] + gt['width'], dt['col'] + dt['width']);
			yy2u = max(gt['row'] + gt['height'], dt['row'] + gt['height']);
			wi = max(float(xx2i - xx1i), 0.0);
			hi = max(float(yy2i - yy1i), 0.0);
			wu = float(xx2u - xx1u);
			hu = float(yy2u - yy1u);
			iou_ratio = wi*hi / (wu*hu);
			# assign dt to gt
			if iou_ratio >= ov_thresh:
				dt_assign[dt_cnt - 1] = gt_cnt - 1;
				detected = True;
	# calculate TP, FP for top K
	for k in range(0,len(sorted_dt_array_symbol)):
		topK_sorted_dt_array_symbol = sorted_dt_array_symbol[:k+1];
		topK_dt_assign = dt_assign[:k+1];
		tp = 0;
		fp = 0;
		# unassigned detection
		fp = len(topK_dt_assign[np.where(topK_dt_assign == -1)]);
		for kk in range(0,len(gt_array_symbol)):
			num_detections = len(topK_dt_assign[np.where(topK_dt_assign == kk)]);
			if num_detections > 0:
				tp = tp + 1;
				fp = fp + num_detections - 1;
		precision = 0.0;
		recall = 0.0;
		if k == len(sorted_dt_array_symbol) - 1:
			total_tp = total_tp + tp;
		if tp + fp > 0:
			precision = float(tp)/float(tp + fp);
		if len(gt_array_symbol) > 0:
			recall = float(tp)/float(len(gt_array_symbol));
		pr_array[k] = precision;
		rc_array[k] = recall;

	if len(pr_array) > 0:
		print "Mean average precision for %s = %f" % (symbol_set[i], np.mean(pr_array));
		mAP = mAP + np.mean(pr_array);
	else:
		print "Mean average precision for %s = %f" % (symbol_set[i], 0);
	# make plots
	pl.clf();
	pl.plot(rc_array, pr_array, label='Precision-Recall curve');
	pl.xlabel('Recall');
	pl.ylabel('Precision');
	pl.ylim([0.0, 1.05]);
	pl.xlim([0.0, 1.0]);
	pl.title('Precision-Recall for %s' % symbol);
	pl.legend(loc="best");
	pl.savefig(symbol+"_PRCurve.png");
	# pl.show();



