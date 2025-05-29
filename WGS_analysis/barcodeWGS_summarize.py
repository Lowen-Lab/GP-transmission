import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

min_read_count = 500

plot_mode = 'stacked_bars' # 'frequency'

project_dir = "./"
input_dir = project_dir + "barcode_count_files/"
filelist = [f for f in os.listdir(input_dir) if f.endswith(".txt")]


barcode_color_dict = {}
infile = open("barcode_color_map.kate.txt","r")
for line in infile:
	line = line.strip().split("\t")
	barcodeID = line[0]
	color = line[1]
	barcode_color_dict[barcodeID] = color
infile.close()


accession_list = []
wgs_barcode_count_dict = {}
wgs_count_dict = {}
for files in filelist:
	accession = files.split('.')[0]
	wgs_count_dict[accession] = 0
	wgs_barcode_count_dict[accession] = {}
	accession_list.append(accession)
	infile = open(input_dir+files,'r')
	first_line = True
	for line in infile:
		if first_line == True:
			first_line = False
		else:
			line = line.strip().split('\t')
			barcode_nt_string = line[0]
			barcode_mismatch_count = int(line[1])
			backbone_mismatch_count = int(line[2])
			if barcode_mismatch_count == 0 and backbone_mismatch_count == 0:
				try:
					color = barcode_color_dict[barcode_nt_string]
				except:
					color = ''
				if color != '':
					wgs_count_dict[accession] += 1
					try:
						wgs_barcode_count_dict[accession][barcode_nt_string] += 1
					except:
						wgs_barcode_count_dict[accession][barcode_nt_string] = 1
	infile.close()

outlines = ''
accessions_to_plot_list = []
for accession in wgs_count_dict:
	count_total = wgs_count_dict[accession]
	outlines += accession+'\t'+str(count_total)+'\n'
	if count_total >= min_read_count:
		accessions_to_plot_list.append(accession)

# outfile = open('wgs_read_count.txt','w')
# outfile.write(outlines)
# outfile.close()
# sys.exit()

full_barcode_list = list(barcode_color_dict.keys())

freq_infile = open("all_barcode_freq.table.txt","r")
first_line = True
freq_dict = {}
sample_to_column_num_dict = {}
barcode_list = []
for line in freq_infile:
	line = line.strip().split("\t")
	if first_line == True:
		sample_list = line
		first_line = False
		for i in range(0,len(sample_list)):
			sampleID = sample_list[i].split('_')[0]
			sample_to_column_num_dict[sampleID] = i
	else:
		barcodeID = line[0]
		barcode_list.append(barcodeID)
		for number in range(1,len(line)):
			freq = float(line[number])
			sampleID = sample_list[number-1].split('_')[0]
			if freq >0 and sampleID in accessions_to_plot_list:
				try:
					freq_dict[sampleID][barcodeID] = freq
				except:
					freq_dict[sampleID] = {}
					freq_dict[sampleID][barcodeID] = freq
freq_infile.close()


num_col = 4
num_row = math.ceil(len(accessions_to_plot_list)/num_col)
col = -1
row = -1

fig, axs = plt.subplots(num_row, num_col, sharex=True, sharey=True)
fig.set_size_inches(num_col*3, num_row*3)


for i in range(0,len(accessions_to_plot_list)):
	sampleID = accessions_to_plot_list[i]
	col +=1
	if i % num_col == 0:
		col = 0
		row += 1
	print(str(row)+' '+str(col))
	count_wgs_total = wgs_count_dict[sampleID]
	local_barcode_list = []
	amp_freq_list = []
	wgs_freq_list = []
	local_color_list = []
	df_key = ['day']
	for b in range(0,len(full_barcode_list)):
		barcodeID = full_barcode_list[b]
		barcode_color = barcode_color_dict[barcodeID]
		try:
			amp_freq = freq_dict[sampleID][barcodeID]
		except:
			amp_freq = 1e-6

		try:
			wgs_count = wgs_barcode_count_dict[sampleID][barcodeID]
			wgs_freq = wgs_count/count_wgs_total
		except:
			wgs_freq = 1e-6

		if wgs_freq >1e-6 or amp_freq>1e-6:
			amp_freq_list.append(amp_freq)
			wgs_freq_list.append(wgs_freq)
			local_color_list.append(barcode_color)
			local_barcode_list.append(barcodeID)
	df_key.extend(local_barcode_list)
	if plot_mode == 'frequency':
		axs[row, col].scatter(amp_freq_list,wgs_freq_list,s=20,linewidths=0)
		axs[row, col].set_xlabel("Amplicon frequency")
		axs[row, col].set_ylabel("WGS frequency")
		axs[row, col].set_xscale('log')
		axs[row, col].set_yscale('log')
		axs[row,col].set_xlim(1e-6,1.5)
		axs[row,col].set_ylim(1e-6,1.5)
	elif plot_mode == 'stacked_bars':
		barcode_array = [amp_freq_list,wgs_freq_list]
		barcode_df = pd.DataFrame(barcode_array,columns=df_key)
		bar_plot = barcode_df.plot(x='day', kind='bar', stacked=True,legend=False,color=local_color_list,ax=axs[row,col],width=0.9)

	title_string = sampleID
	axs[row, col].set_title(title_string)
# plt.tight_layout()
if plot_mode == 'frequency':
	plt.savefig("wgs_to_amplicon_freq.pdf")
elif plot_mode == 'stacked_bars':
	plt.savefig("wgs_ampl_barplots.pdf")
