import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
import pandas
from scipy.stats import gaussian_kde


project_dir = "./"
input_read_dir =project_dir+"raw_data/"
ref_dir = project_dir+"refs/"
trim_dir = project_dir+"trim/"
map_folder = 'rep_map'
map_dir = project_dir+map_folder+"/"
qual_dir = map_dir+"qual/"
sam_dir = map_dir+"sam/"
bam_dir = map_dir+"bam/"
sub_dir = map_dir+"substitutions/"
temp_dir = map_dir+"temp/"

output_suffix = "wgs"

all_segment_list = ['IAV_PB2','IAV_PB1','IAV_PA','IAV_HA','IAV_NP','IAV_NA','IAV_M','IAV_NS']

min_frequency = 0.02
min_cov = 1000
min_minor_qual = 25


num_plots = 6
x = {}
y = {}
for r in range(0,num_plots+1):
	x[r] = {}
	y[r] = {}
	for i in range(0,len(all_segment_list)):
		segment = all_segment_list[i]
		x[r][segment] = []
		y[r][segment] = []

infile = open("all_sites.base_info."+output_suffix+".txt","r")
first_line = True
for line in infile:
	if first_line == True:
		first_line = False
	else:
		line = line.strip().split("\t")
		sampleID = line[0]
		segment = line[1]
		
		loc = int(line[2])
		cov = int(line[3])
		if cov >= min_cov:
			major_prop = round(float(line[6]),1)
			major_qual = round(float(line[7]),1)
			major_mapq = round(float(line[8]),1)
			major_loc = round(float(line[9]),1)
			major_mismatch = round(float(line[10]),1)
			major_indel = round(float(line[11]),1)

			minor_prop = round(float(line[12]),1)
			minor_qual = round(float(line[13]),1)
			minor_mapq = round(float(line[14]),1)
			minor_loc = round(float(line[15]),1)
			minor_mismatch = round(float(line[16]),1)
			minor_indel = round(float(line[17]),1)

			if major_mapq >=30:

				xval = major_mapq
				yval = major_loc
				r = 0
				x[r][segment].append(xval)
				y[r][segment].append(yval)

				xval = loc
				yval = major_qual
				r = 1
				x[r][segment].append(xval)
				y[r][segment].append(yval)

				xval = loc
				yval = major_mismatch
				r = 2
				x[r][segment].append(xval)
				y[r][segment].append(yval)

				xval = loc
				yval = major_indel
				r = 3
				x[r][segment].append(xval)
				y[r][segment].append(yval)

				xval = loc
				yval = major_mapq
				r = 4
				x[r][segment].append(xval)
				y[r][segment].append(yval)

				if minor_prop >= min_frequency and minor_qual >= min_minor_qual:
					xval = minor_mapq
					yval = minor_loc
					r = 5
					x[r][segment].append(xval)
					y[r][segment].append(yval)

infile.close()


for r in range(0,num_plots+1):
	lineout = ''
	for i in range(0,len(all_segment_list)):
		segment = all_segment_list[i]
		lineout += "\t"+str(len(x[r][segment]))
	print(lineout)


fig, axs = plt.subplots(num_plots,len(all_segment_list), sharex=False, sharey=False)
fig.set_size_inches(len(all_segment_list)*4.0,num_plots*3.5)

axis_labels = {0:("mapq","read_loc",(30,50),(0,40)),1:("loc","qual",(0,2400),(30,50)),2:("loc","mismatch",(0,2400),(0,5)),3:("loc","indel",(0,2400),(0,5)),4:("loc","mapq",(0,2400),(30,50)),5:("minor_mapq","minor_read_loc",(20,45),(0,50))}

for r in range(0,num_plots):
	row = r
	for i in range(0,len(all_segment_list)):
		segment = all_segment_list[i]
		if len(x[r][segment]) >0:
			if len(x[r][segment]) > 100000:
				x[r][segment] = x[r][segment][0:100000]
				y[r][segment] = y[r][segment][0:100000]
			print(str(row)+"\t"+str(i))
			xy = np.vstack([x[r][segment],y[r][segment]])
			try:
				z = gaussian_kde(xy)(xy)
				axs[row, i].scatter(x[r][segment],y[r][segment], c=z, s=3)
			except:
				axs[row, i].scatter(x[r][segment],y[r][segment], s=3)

			axs[row, i].set_xlabel(axis_labels[r][0])
			axs[row, i].set_ylabel(axis_labels[r][1])
			if row == 0:
				axs[row, i].set_title(segment)

			axs[row,i].set_xlim([axis_labels[r][2][0], axis_labels[r][2][1]])
			axs[row,i].set_ylim([axis_labels[r][3][0], axis_labels[r][3][1]])

plt.tight_layout()
plt.savefig("qual_summary.major_allele."+str(output_suffix)+".png")
