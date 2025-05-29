import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import math
import pandas as pd

window_size = 5
step_size = 2
nameout_string = "wgs"
project_dir = "./"
ref_dir = project_dir+"refs/"

all_segment_list = ['IAV_PB2','IAV_PB1','IAV_PA','IAV_HA','IAV_NP','IAV_NA','IAV_M','IAV_NS']

def mean_confidence_interval(data, confidence=0.99):
	x = []
	for d in data:
		if d > 0 and isNumber(d) == True:
			x.append(d)
	a = 1.0 * np.array(x)
	m = np.nanmean(data)
	ll = np.nanquantile(data,0.1)
	l = np.nanquantile(data,0.25)
	u = np.nanquantile(data,0.75)
	uu = np.nanquantile(data,0.9)

	outliers = []
	# for i in range(0,len(data)):
	# 	if data[i] > uu or data[i] < ll:
	# 		outliers.append(data[i])
	return m, l, u, outliers, ll, uu

def isNumber(s):
	try:
		float(s)
	except ValueError:
		return False
	return True

def count_in_list(list_in,value=np.nan):
	counter = 0
	for item in list_in:
		if item == value:
			counter += 1
		else:
			try:
				if item is value:
					counter += 1
			except:
				pass
	return counter

consensus_seq_dict = {}
max_length_dict = {}


allcov_infile = open(project_dir+"all_cov."+nameout_string+".txt","r")
first_line = True
cov_dict = {}
max_cov_dict = {}
for line in allcov_infile:
	line = line.strip().split("\t")
	if first_line == True:
		sample_header = line
		first_line = False
	else:
		segment = line[0]
		try:
			loc = int(line[1])
		except:
			print(line)
			sys.exit()
		cov_list = []
		try:
			if loc > max_length_dict[segment]:
				max_length_dict[segment] = loc
		except:
			max_length_dict[segment] = loc
		for col in range(2,len(line)):
			sampleID = sample_header[col]
			depth = int(float(line[col]))
			if depth > 0:
				try:
					cov_dict[segment][sampleID][loc] = depth
				except:
					try:
						cov_dict[segment][sampleID] = {}
						cov_dict[segment][sampleID][loc] = depth
					except:
						cov_dict[segment] = {}
						cov_dict[segment][sampleID] = {}
						cov_dict[segment][sampleID][loc] = depth
allcov_infile.close()

cov_list_dict = {}
cov_bin_dict = {}
for segment in cov_dict:
	for sampleID in cov_dict[segment]:
		cov_list = []
		for loc in range(0,max_length_dict[segment]):
			try:
				cov = cov_dict[segment][sampleID][loc]
			except:
				cov = np.nan
			cov_list.append(cov)
		
		try:
			cov_list_dict[segment][sampleID] = cov_list
		except:
			cov_list_dict[segment] = {}
			cov_list_dict[segment][sampleID] = cov_list

avg_cov_dict = {}
seg_loc_dict = {}
stat_dict = {}
for i in range(0,len(all_segment_list)):
	segment = all_segment_list[i]
	seg_loc_dict[segment] = []
	print(segment)
	for loc in range(0,max_length_dict[segment]-int(window_size/2),step_size):
		low = int(loc-int(window_size-1)/2)
		high = int(loc+int(window_size-1)/2)+1
		seg_loc_dict[segment].append(loc)
		for sampleID in cov_dict[segment]:
			try:
				cov_list = cov_list_dict[segment][sampleID][low:high]
			except:
				cov_list = []
			if cov_list != []:
				if count_in_list(cov_list) == 0:
					avg_cov = np.nanmean(cov_list)
					try:
						avg_cov_dict[segment][sampleID][loc] = avg_cov
					except:
						try:
							avg_cov_dict[segment][sampleID] = {}
							avg_cov_dict[segment][sampleID][loc] = avg_cov
						except:
							avg_cov_dict[segment] = {}
							avg_cov_dict[segment][sampleID] = {}
							avg_cov_dict[segment][sampleID][loc] = avg_cov

					try:
						stat_dict[segment][loc].append(avg_cov)
					except:
						try:
							stat_dict[segment][loc] = []
							stat_dict[segment][loc].append(avg_cov)
						except:
							try:
								stat_dict[segment] = {}
								stat_dict[segment][loc] = []
								stat_dict[segment][loc].append(avg_cov)
							except:
								stat_dict = {}
								stat_dict[segment] = {}
								stat_dict[segment][loc] = []
								stat_dict[segment][loc].append(avg_cov)

lineout = ''
x = {}
mean = {}
u95 = {}
l95 = {}
u99 = {}
l99 = {}
outlier_dict = {}
outlier_loc_dict = {}
for i in range(0,len(all_segment_list)):
	segment = all_segment_list[i]
	try:
		num = len(stat_dict[segment])
	except:
		num = 0
	lineout += "\t"+str(num)
	if num >0:
		x[segment] = []
		mean[segment] = []
		u95[segment] = []
		l95[segment] = []
		u99[segment] = []
		l99[segment] = []
		outlier_dict[segment] = []
		outlier_loc_dict[segment] = []
		for loc in range(0,max_length_dict[segment]-int(window_size/2),step_size):
			try:
				stat_list = stat_dict[segment][loc]
				m,l,u,outliers, ll, uu = mean_confidence_interval(stat_list)
				x[segment].append(loc)
				mean[segment].append(m)
				u95[segment].append(u)
				l95[segment].append(l)
				u99[segment].append(uu)
				l99[segment].append(ll)
			except:
				pass
print(lineout)


fig, axs = plt.subplots(1,len(all_segment_list), sharex=True, sharey=True)
fig.set_size_inches(len(all_segment_list)*3,3)

for i in range(0,len(all_segment_list)):
	segment = all_segment_list[i]
	try:
		num = len(stat_dict[segment])
	except:
		num = 0
	if num >0:
		loc_list = x[segment]
		mean_list = mean[segment]
		u95_list = u95[segment]
		l95_list = l95[segment]
		u99_list = u99[segment]
		l99_list = l99[segment]
		seg_outname = segment
		
		axs[i].fill_between(loc_list, u95_list, l95_list, color = 'dodgerblue', alpha = 1.0)
		axs[i].fill_between(loc_list, u95_list, u99_list, color = 'deepskyblue', alpha = 1.0)
		axs[i].fill_between(loc_list, l95_list, l99_list, color = 'deepskyblue', alpha = 1.0)
		axs[i].plot(loc_list,mean_list,color="red", linewidth=0.5)
		axs[i].set_title(seg_outname)
		axs[i].set_yscale("log")
		axs[i].set_ylim([1e0, 1e6])
plt.tight_layout()
plt.savefig("cov_summary."+nameout_string+".pdf")
