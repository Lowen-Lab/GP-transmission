import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

project_dir = "./"
number_of_files_split_figures = 1
nameout_string = "titer_scaled" #"not_scaled"
fig_name_prefix = "barcode_summary"
min_titer = 1.5
min_freq_to_plot = 0.

sample_info_dict = {}
gp_info_dict = {}
pair_dict = {}
titer_dict = {}
pair_list = []
max_titer = 0
filename = "sample_key.txt"
infile = open(project_dir+filename,"r")
first_line = True
for line in infile:
	if first_line == True:
		first_line = False
	else:
		line = line.strip().split("\t")
		GPID = line[7]
		dayID = line[5]
		day = int(line[4])
		sampleID = line[1]
		pairID = line[9]
		GP_type = line[10]
		route = line[11]
		exp = line[2]
		titer = float(line[12])
		titer_dict[sampleID] = titer
		if GP_type != "s": 
			if titer > max_titer:
				max_titer = titer
			try:
				pair_dict[pairID][GP_type] = GPID
			except:
				pair_dict[pairID] = {}
				pair_dict[pairID][GP_type] = GPID
			pair_list.append(pairID)

			tup = (GPID,pairID,GP_type,route)
			gp_info_dict[GPID] = tup

			tup = (GPID,dayID,sampleID,pairID,GP_type,route,exp,titer)
			try:
				sample_info_dict[GPID][day] = tup
			except:
				sample_info_dict[GPID] = {}
				sample_info_dict[GPID][day] = tup
infile.close()
pair_list = sorted(list(set(pair_list)))


barcode_color_dict = {}
infile = open("barcode_color_map.txt","r")
for line in infile:
	line = line.strip().split("\t")
	barcode = line[0]
	color = line[1]
	barcode_color_dict[barcode] = color
infile.close()


freq_infile = open("./barcode_info/all_barcode_freq.table.txt","r")
first_line = True
freq_dict = {}
color_list = []
sample_to_column_num_dict = {}
barcode_list = []
for line in freq_infile:
	line = line.strip().split("\t")
	if first_line == True:
		sample_list = line
		first_line = False
		for i in range(0,len(sample_list)):
			sampleID = sample_list[i]
			sample_to_column_num_dict[sampleID] = i
	else:
		barcodeID = line[0]
		barcode_color = barcode_color_dict[barcodeID]
		color_list.append(barcode_color)
		barcode_list.append(barcodeID)
		for number in range(1,len(line)):
			freq = float(line[number])
			if freq <= min_freq_to_plot:
				freq = np.nan
			sampleID = sample_list[number-1]
			try:
				freq_dict[sampleID].append(freq)
			except:
				freq_dict[sampleID] = []
				freq_dict[sampleID].append(freq)
freq_infile.close()


### make arrays for stacked bar plots
adjust_freq_dict = {}
for sampleID in freq_dict:
	titer = titer_dict[sampleID]
	freq_array = freq_dict[sampleID]
	adjust_freq_dict[sampleID] = []
	for f in range(0,len(freq_array)):
		if nameout_string == "not_scaled":
			adjust_freq_dict[sampleID].append(freq_array[f])
		else:
			adjust_freq_dict[sampleID].append(freq_array[f]*titer)


print("number of infection pairs to plot: "+str(len(pair_list)))
num_col = 6
num_row = int((len(pair_list)/number_of_files_split_figures))
col = 0 #column
row = -1 #row
fig_num = -1
used_pigIDs = []
print("\tpairs per plot: "+str(num_row))

fig, axs = plt.subplots(int(num_row/3), num_col, sharex=False, sharey=False)
fig.set_size_inches(num_col*2.0, num_row*2.5/3)

for r in range(0,len(pair_list)):
	pairID = pair_list[r]
	GP_i = pair_dict[pairID]["i"]
	GP_e = pair_dict[pairID]["e"]
	print(pairID)
	try:
		gp_i_info = gp_info_dict[GP_i] #(GPID,pairID,GP_type,route)
	except:
		gp_i_info = ()
	try:
		gp_e_info = gp_info_dict[GP_e]
	except:
		gp_e_info = ()

	### Save current figure and start new figure if number of rows%rows-per-fig == 0
	if r % num_row == 0 and r > 0:
		row = 0
		plt.tight_layout()
		fig_num += 1
		print("Saving plot")
		plt.savefig(fig_name_prefix+"."+nameout_string+"."+str(fig_num)+".pdf")
		fig, axs = plt.subplots(num_row, num_col, sharex=False, sharey=False)
		fig.set_size_inches(num_col*3.5, num_row*3)
	else:
		if r % 3 == 0:
			row += 1
			col = -1

	### Make scaled bar plots with barcodes scaled to sample titers
	gp_i_barcode_dict = []
	gp_e_barcode_dict = []
	gp_i_max_mismatch = []
	gp_e_max_mismatch = []
	df_key = ['day']
	df_key.extend(barcode_list)
	for day in range(1,9):
		try:
			i_sampleID = sample_info_dict[GP_i][day][2]
			adjust_freq_dict[i_sampleID]
		except:
			i_sampleID = ''
		if i_sampleID != '':
			i_titer = titer_dict[i_sampleID]
			if i_titer >= min_titer:
				i_sampleID = ''
		if i_sampleID != '':
			i_titer_freq_array = [day]
			for i in range(0,len(adjust_freq_dict[i_sampleID])):
				i_titer_freq_array.append(adjust_freq_dict[i_sampleID][i])
			gp_i_barcode_dict.append(i_titer_freq_array)

		try:
			e_sampleID = sample_info_dict[GP_e][day][2]
			adjust_freq_dict[e_sampleID]
		except:
			e_sampleID = ''
		if e_sampleID != '':
			e_titer = titer_dict[e_sampleID]
			if e_titer >= min_titer:
				e_sampleID = ''
		if e_sampleID != '':
			e_titer_freq_array = [day]
			for i in range(0,len(adjust_freq_dict[e_sampleID])):
				e_titer_freq_array.append(adjust_freq_dict[e_sampleID][i])
			gp_e_barcode_dict.append(e_titer_freq_array)
	
	print("\t"+GP_i)
	col += 1
	if gp_i_barcode_dict != []:
		gp_i_barcode_df = pd.DataFrame(gp_i_barcode_dict,columns=df_key)
		title_string = GP_i+" - i - "+gp_i_info[1]+" - "+gp_i_info[3]
		if nameout_string == "not_scaled":
			bar_plot_i = gp_i_barcode_df.plot(x='day', kind='bar', stacked=True,legend=False,color=color_list,ax=axs[row,col],width=0.9,title=title_string, ylim=(0, 1.05),ylabel="Frequency")# , ylim=(0, max_titer*1.05),ylabel="Titer (pfu/mL)")
		else:
			bar_plot_i = gp_i_barcode_df.plot(x='day', kind='bar', stacked=True,legend=False,color=color_list,ax=axs[row,col],width=0.9,title=title_string, ylim=(0, max_titer*1.05),ylabel="Titer (pfu/mL)")#, ylim=(0, 1.05),ylabel="Frequency")
			bar_plot_i.axhline(y=min_titer, color = 'red', zorder=-1)
		bar_plot_i.set_xticks(range(1,10))

	print("\t"+GP_e)
	col += 1
	if gp_e_barcode_dict != []:
		gp_e_barcode_df = pd.DataFrame(gp_e_barcode_dict,columns=df_key)
		title_string = GP_e+" - e - "+gp_i_info[1]+" - "+gp_i_info[3]
		if nameout_string == "not_scaled":
			bar_plot_e = gp_e_barcode_df.plot(x='day', kind='bar', stacked=True,legend=False,color=color_list,ax=axs[row,col],width=0.9,title=title_string, ylim=(0, 1.05),ylabel="Frequency")#, ylim=(0, max_titer*1.05),ylabel="Titer (pfu/mL)")
		else:
			bar_plot_e = gp_e_barcode_df.plot(x='day', kind='bar', stacked=True,legend=False,color=color_list,ax=axs[row,col],width=0.9,title=title_string, ylim=(0, max_titer*1.05),ylabel="Titer (pfu/mL)")#, ylim=(0, 1.05),ylabel="Frequency")
			bar_plot_e.axhline(y=min_titer, color = 'red', zorder=-1)
		bar_plot_e.set_xticks(range(1,10))

fig_num += 1
plt.tight_layout()
print("Saving plot")
plt.savefig(fig_name_prefix+"."+nameout_string+"."+str(fig_num)+".pdf")
# plt.show()
