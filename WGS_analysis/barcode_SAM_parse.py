import os
import sys
import numpy as np
import multiprocessing
from joblib import Parallel, delayed

dedupe_reads = True
skip_qual = False

parallel_max_cpu = 30

###################################   Setup   ###################################

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
barcode_dir = map_dir+"barcodes/"
command_prefix = "bash "
ref_set_list = ['IAV_NA.fa']
focal_region_start = 483
focal_region_end = 531
barcode_loc_dict = {0:(483,'T,C'),1:(486,'A,G'),2:(490,'C,T'),3:(504,'A,G'),4:(507,'A,G'),5:(513,'T,C'),6:(519,'T,C'),7:(522,'T,C'),8:(523,'T,C'),9:(525,'A,G'),10:(528,'A,G'),11:(531,'T,C')}
backbone_loc_dict = {427:'G',428:'G',429:'A',430:'A',431:'C',432:'A',433:'A',434:'C',435:'A',436:'C',437:'T',438:'A',439:'A',440:'A',441:'C',442:'A',443:'A',444:'C',445:'A',446:'G',447:'G',448:'C',449:'A',450:'T',451:'T',452:'C',453:'A',454:'A',455:'A',456:'T',457:'G',458:'A',459:'C',460:'A',461:'C',462:'A',463:'G',464:'T',465:'A',466:'C',467:'A',468:'T',469:'G',470:'A',471:'T',472:'A',473:'G',474:'G',475:'A',476:'C',477:'C',478:'C',479:'C',480:'T',481:'T',482:'A',484:'C',485:'G',487:'A',488:'C',489:'C',491:'T',492:'A',493:'T',494:'T',495:'G',496:'A',497:'T',498:'G',499:'A',500:'A',501:'T',502:'G',503:'A',505:'T',506:'T',508:'G',509:'G',510:'T',511:'G',512:'T',514:'C',515:'C',516:'A',517:'T',518:'T',520:'C',521:'A',524:'T',526:'G',527:'G',529:'A',530:'C',532:'A',533:'A',534:'G',535:'C',536:'A',537:'A',538:'G',539:'T',540:'G',541:'T',542:'G',543:'T',544:'A',545:'T',546:'A',547:'G',548:'C',549:'A',550:'T',551:'G',552:'G',553:'T',554:'C',555:'C',556:'A',557:'G',558:'C',559:'T',560:'C',561:'A',562:'A',563:'G',564:'T',565:'T',566:'G',567:'T',568:'C',569:'A',570:'C',571:'G',572:'A',573:'T',574:'G',575:'G',576:'A',577:'A',578:'A',579:'A',580:'G',581:'C'}

###################################   Setup   ###################################
#Set up parallel processing
if parallel_max_cpu == 0:
	num_cores = multiprocessing.cpu_count()
elif parallel_max_cpu > 0:
	num_cores = min(parallel_max_cpu,multiprocessing.cpu_count())
elif parallel_max_cpu < 0:
	num_cores = multiprocessing.cpu_count()+parallel_max_cpu
else:
	num_cores = 1

#Make output directories if they don't exist yet
if not os.path.exists(barcode_dir):
	os.makedirs(barcode_dir)
if not os.path.exists(sub_dir):
	os.makedirs(sub_dir)
if not os.path.exists(temp_dir):
	os.makedirs(temp_dir)
if skip_qual == False:
	if not os.path.exists(qual_dir):
		os.makedirs(qual_dir)

############################################ FUNCTIONS ############################################
def run_command(command,mode="quiet"):
	run_status = False
	if mode == "verbose":
		return_code = os.system(command)
	elif mode == "quiet":
		return_code = os.system(command+" >/dev/null 2>&1")
	if return_code == 0:
		run_status = True
	return run_status

def shannon_alpha(freq_array):
	H = 0
	S_obs = 0
	for freq in freq_array:
		if freq >0.0:
			H_i = -1*freq*np.log(freq)
			H += H_i
			S_obs += 1
	H_max = np.log(S_obs)
	E = H/H_max
	H = round(H,3)
	E = float(round_to_n_sig_figs(E,3))
	return H,E

def simpson_alpha(freq_array):
	p_sum = 0
	S = 0
	for freq in freq_array:
		if freq >0.0:
			S += 1
			p_sum += freq**2
	H = p_sum
	return round_to_n_sig_figs(H,3)

def bray_curtis_dissimilarity_log(array1,array2):
	C,S1,S2 = 0,0,0
	for i in range(0,len(array1)):
		val1 = array1[i]
		val2 = array2[i]
		if val1 > 0.0 and val2 > 0.0:
			minval = min(val1,val2)
			C += 1/(-1*np.log(min(val1,val2)))
		if val1 > 0.0:
			S1 += 1/(-1*np.log(val1))
		if val2 > 0.0:
			S2 += 1/(-1*np.log(val2))
	B = 1-((2*C)/(S1+S2))
	return round_to_n_sig_figs(B,3)

def bray_curtis_dissimilarity(array1,array2):
	C,S1,S2 = 0,0,0
	for i in range(0,len(array1)):
		val1 = array1[i]
		val2 = array2[i]
		if val1 > 0.0 and val2 > 0.0:
			C += min(val1,val2)
		if val1 > 0.0:
			S1 += val1
		if val2 > 0.0:
			S2 += val2
	B = 1-((2*C)/(S1+S2))
	return round_to_n_sig_figs(B,3)

def chao_1_richness(count_array):
	S_obs = 0
	S_single = 0
	S_double = 0
	for count in count_array:
		if count >0:
			S_obs += 1
			if count == 1:
				S_single += 1
			elif count == 2:
				S_double += 1
	S_unobs = (S_single*(S_single-1))/(2*S_double+1)
	R = S_obs + S_unobs
	return round(R,1)

def is_number(n):
	try:
		float(n)
	except ValueError:
		return False
	return True

def round_to_n_sig_figs(val,num_sig_figs):
	if is_number(val):
		val = float(val)
		if val == 0.0:
			return '0.0'
		num_sig_figs = int(num_sig_figs)
		if num_sig_figs == 0:
			num_sig_figs = 1
		sci_val = "{:.10e}".format(val)
		split_sci_val = sci_val.split("e")
		if len(split_sci_val) == 2:
			rounded_base_number = round(float(split_sci_val[0]),num_sig_figs-1)
			exponent = int(split_sci_val[1])
			if exponent == 0:
				val_out = str(rounded_base_number) + ((num_sig_figs)-1)*'0'
			elif exponent < 0:
				exponent*=-1
				val_out = '0.' + (exponent-1)*'0' + str(rounded_base_number).replace(".","")
				val_out = str(float(val_out))
			elif exponent > 0:
				val_out = str(rounded_base_number) +'e'+ (str(exponent))
			return val_out
		else:
			return val
	else:
		return val


def check_if_reads_mapped(sampleID,ref_set,sam_dir,read_count_dict):
	reads_mapped = True
	refs = ref_set.split("-")
	primary_ref = ref_set.split("-")[0].split(".f")[0]
	for ref_filename in refs:
		if os.path.isfile(sam_dir+primary_ref+"/"+sampleID+"."+ref_filename.split(".f")[0]+'.sam') == False:
			if read_count_dict[ref_filename.split(".f")[0]] <= 100:
				reads_mapped = False
	return reads_mapped


def check_if_processing_finished(sampleID,ref_set,sub_dir):
	sam_processed = True
	refs = ref_set.split("-")
	primary_refID = ref_set.split("-")[0].split(".f")[0]
	out_substitutions_filename = sub_dir+primary_refID+"/"+sampleID+"."+primary_refID+".substitions.txt"
	if os.path.isfile(out_substitutions_filename) == False:
		sam_processed = False
	return sam_processed


def load_read_counts(sampleID,trim_dir):
	read_count_dict = {}
	for ref_set in ref_set_list:
		refs = ref_set.split("-")
		primary_ref = ref_set.split("-")[0].split(".f")[0]
		for ref_filename in refs:
			segment = ref_filename.split(".f")[0]
			read_count_file = open(trim_dir+primary_ref+"/"+sampleID+'.'+segment+'.read_count.txt','r')
			for line in read_count_file:
				line = line.strip()
				read_count_dict[segment] = int(line)
			read_count_file.close()
	return read_count_dict


def dedupe_samfile(sam_filename_in,temp_sam_filename,sam_filename_out):
	sam_infile = open(sam_filename_in,"r")
	read_map_strings = {}
	for sam_line in sam_infile:
		if len(sam_line) > 0:
			if sam_line[0] != "@":
				sam_line = sam_line.strip().split("\t")
				read_seq = sam_line[9]
				if read_seq !="*":
					readID = sam_line[0].split(" ")[0]
					map_site = sam_line[3]
					map_qual = sam_line[4]
					cigar = sam_line[5]
					map_string = map_site+"_"+map_qual+"_"+cigar
					try:
						read_map_strings[readID] += " "+map_string
					except:
						read_map_strings[readID] = " "+map_string
	sam_infile.close()
	map_string_dedupe = {}
	for readID in read_map_strings:
		map_string = read_map_strings[readID]
		if len(map_string.split(" ")[0]) == 2:
			map_string_R = map_string.split(" ")[0]+" "+map_string.split(" ")[1]
			try:
				one_readID = map_string_dedupe[map_string]
			except:
				try:
					one_readID = map_string_dedupe[map_string_R]
				except:
					map_string_dedupe[map_string] = readID
		else:
			try:
				one_readID = map_string_dedupe[map_string]
			except:
				map_string_dedupe[map_string] = readID
	del read_map_strings
	pass_readID_dict = {}
	for map_string in map_string_dedupe:
		readID = map_string_dedupe[map_string]
		pass_readID_dict[readID] = ''
	del map_string_dedupe
	dedupe_sam_outfile = open(temp_sam_filename,"w")
	sam_infile = open(sam_filename_in,"r")
	for sam_line in sam_infile:
		if len(sam_line) > 0:
			sam_line = sam_line.strip()#.split("\t")
			if sam_line[0] != "@":
				split_line = sam_line.split("\t")
				readID = split_line[0].split(" ")[0]
				read_seq = split_line[9]
				if read_seq !="*":
					try:
						pass_readID_dict[readID]
						dedupe_sam_outfile.write(sam_line+"\n")
					except:
						pass
			else:
				dedupe_sam_outfile.write(sam_line+"\n")
	sam_infile.close()
	dedupe_sam_outfile.close()
	os.rename(temp_sam_filename,sam_filename_out)
	return len(pass_readID_dict)


def SAM_parse_barcode(map_position,cigar_string,quality_string,read_seq,mapping_qual,focal_region_start,focal_region_end,barcode_loc_dict,backbone_loc_dict): #mapping position is indexed to 1, not zero
	cigar_dict = {'M':'imperfect_match','D':'insertion_in_ref','I':'insertion_in_read','S':'soft_clip_reference','H':'soft_clip_read','=':'perfect_match','X':'mismatch'}
	
	read_key = ''
	ref_site_map = []
	match_count_dict = {}
	
	cur_position = int(map_position)-2
	read_start_site = cur_position
	cur_num = ''
	cur_operator = ''
	all_valid_operators = True
	for i in range(0,len(cigar_string)):
		character = cigar_string[i]
		try:
			cigar_dict[character]
			cur_operator = character

			cur_num = int(cur_num)
			try:
				match_count_dict[cur_operator] += cur_num
			except:
				match_count_dict[cur_operator] = cur_num
			for num in range(0,cur_num):
				if cur_operator == "D":
					cur_position += 1
				elif cur_operator == "I":
					cur_position += 0
					ref_site_map.append('')
					read_key += cur_operator
				elif cur_operator == "S":
					cur_position += 0
					ref_site_map.append('')
					read_key += cur_operator
				elif cur_operator == "H":
					all_valid_operators = False
				elif cur_operator == "=": #perfect match
					cur_position += 1
					ref_site_map.append(cur_position)
					read_key += cur_operator
				elif cur_operator == "M": #imperfect match
					cur_position += 1
					ref_site_map.append(cur_position)
					read_key += cur_operator
				elif cur_operator == "X": #mismatch
					cur_position += 1
					ref_site_map.append(cur_position)
					read_key += cur_operator
				else:
					all_valid_operators = False
					sys.exit("Invalid CIGAR operator encountered. Exiting.")
			cur_num = ''
			cur_operator = ''
		except:
			cur_num += character
	base_dict_fun = {}
	read_end_site = cur_position
	
	if read_end_site >= focal_region_end and read_start_site <= focal_region_start:
		for num in range(0,len(read_key)):
			operator = read_key[num]
			ref_loc = ref_site_map[num]
			if ref_loc != '':
				ref_loc = int(ref_loc)
				try:
					qual_char = quality_string[num]
				except:
					sys.exit("Index ")
				qual_score = ord(qual_char)-33
				read_nt = read_seq[num]
				read_length = len(read_seq)
				nt_prop_location = round(num/read_length,2)
				min_end_dist = min(read_length-num,read_length-num)
				tup = (read_nt,qual_score,mapping_qual,read_length,num,min_end_dist,nt_prop_location)
				base_dict_fun[ref_loc] = tup

		barcode_base_string = ''
		barcode_match_string = ''
		barcode_qual_string = ''
		barcode_mismatches = 0
		for b in range(0,len(barcode_loc_dict)):
			barcode_tup = barcode_loc_dict[b]
			ref_loc = barcode_tup[0]
			barcode_bases = barcode_tup[1].split(",")

			try:
				read_nt = base_dict_fun[ref_loc][0]
				qual_score = base_dict_fun[ref_loc][1]
				mapping_qual = base_dict_fun[ref_loc][2]
				read_length = base_dict_fun[ref_loc][3]
				num = base_dict_fun[ref_loc][4]
				min_end_dist = base_dict_fun[ref_loc][5]
				nt_prop_location = base_dict_fun[ref_loc][6]
			except:
				return {},{},'','','','','','','',''

			barcode_base_string += read_nt
			qual_pass = True
			base_pass = True
			if qual_score < 25:
				qual_pass = False

			if qual_pass == False:
				barcode_qual_string += 'X'
			else:
				barcode_qual_string += 'O'

			if read_nt not in barcode_bases:
				barcode_match_string += "X"
				base_pass = False
			else:
				barcode_match_string += "O"
			

		backbone_base_string = ''
		backbone_match_string = ''
		backbone_qual_string = ''
		backbone_mismatches = 0
		for ref_loc in range(427,582):
			try:
				ref_base = backbone_loc_dict[ref_loc]
			except:
				ref_base = '' #if it can't find the site in the non-barcode location dict, the site is a barcode location
				backbone_base_string += "."
				backbone_match_string += "."
				backbone_qual_string += "."
			if ref_base != '':
				try:
					read_nt = base_dict_fun[ref_loc][0]
					qual_score = base_dict_fun[ref_loc][1]
					mapping_qual = base_dict_fun[ref_loc][2]
					read_length = base_dict_fun[ref_loc][3]
					num = base_dict_fun[ref_loc][4]
					min_end_dist = base_dict_fun[ref_loc][5]
					nt_prop_location = base_dict_fun[ref_loc][6]
					
					backbone_base_string += read_nt
					qual_pass = True
					base_pass = True

					if qual_score < 25:
						if min_end_dist >= 10:
							qual_pass = False
					if qual_pass == False:
						backbone_qual_string += 'X'
					else:
						backbone_qual_string += 'O'
					

					if read_nt != ref_base:
						backbone_match_string += "X"
						base_pass = False
					else:
						backbone_match_string += "O"
				except:
					backbone_base_string += "-"
					backbone_qual_string += "-"
					backbone_match_string += "-"
		return base_dict_fun,match_count_dict,barcode_base_string,barcode_match_string,barcode_qual_string,barcode_mismatches,backbone_base_string,backbone_match_string,backbone_qual_string,backbone_mismatches
	else:
		return {},{},'','','','','','','',''

###########################################################

def main_pipeline(sample_ref_tup,ref_set_list,input_read_dir,ref_dir,barcode_dir,trim_dir,map_dir,qual_dir,sam_dir,bam_dir,temp_dir,sub_dir,command_prefix,dedupe_reads,focal_region_start,focal_region_end,barcode_loc_dict,backbone_loc_dict):
	sampleID = sample_ref_tup[0]
	ref_set = sample_ref_tup[1]
	refs = ref_set.split("-")
	primary_refID = refs[0].split(".f")[0]

	read_count_dict = load_read_counts(sampleID,trim_dir)
	reads_mapped = check_if_reads_mapped(sampleID,ref_set,sam_dir,read_count_dict)
	if reads_mapped == False:
		return None
	
	read_map_info_dict = {}
	map_qual_dict = {}
	indel_dict = {}
	mismatch_dict = {}
	for ref_filename in refs:
		refID = ref_filename.split(".f")[0]
		sam_filename = sam_dir+primary_refID+"/"+sampleID+"."+ref_filename.split(".f")[0]+".refine.sam"
		cont = True
		if dedupe_reads  == False:
			sam_filename_to_open = sam_filename
			if os.path.isfile(sam_filename) == False:
				cont = False
		elif dedupe_reads == True:
			if os.path.isfile(sam_filename) == True:
				temp_sam_filename = temp_dir+sampleID+"."+ref_filename.split(".f")[0]+".dedupe.sam"
				dedupe_samfilename = sam_dir+ref_filename.split(".f")[0]+"/"+sampleID+"."+ref_filename.split(".f")[0]+".dedupe.sam"
				num_reads_after_dedupe = dedupe_samfile(sam_filename,temp_sam_filename,dedupe_samfilename)
				sam_filename_to_open = dedupe_samfilename
			else:
				cont = False
		if cont == True:
			print(sample_ref_tup)
			sam_infile = open(sam_filename_to_open,"r")
			temp_barcode_filename = temp_dir+sampleID+"."+ref_filename.split(".f")[0]+".barcodes.txt"
			barcode_filename_out = barcode_dir+primary_refID+"/"+sampleID+"."+ref_filename.split(".f")[0]+".barcodes.txt"
			barcode_outfile = open(temp_barcode_filename,"w")
			barcode_outfile.write('barcode_base_string\tbarcode_mismatch_count\tbackbone_mismatch_count\tbarcode_qual_string\tbarcode_match_string\tbackbone_base_string\tbackbone_qual_string\tbackbone_match_string\tinsertion_count\tdeletion_count\n')
			barcode_region_reads_found = 0
			barcode_reads_found = 0
			barcode_count_dict = {}
			for sam_line in sam_infile:
				if len(sam_line) > 0:
					if sam_line[0] != "@":
						sam_line = sam_line.strip().split("\t")
						map_qual = int(sam_line[4])
						read_seq = sam_line[9]
						if read_seq !="*":
							readID = sam_line[0].split(" ")[0]
							# flag = sam_line[1]
							# contig = sam_line[2]
							map_site = int(sam_line[3])
							cigar = sam_line[5] #{'M':'imperfect_match','D':'insertion_in_ref','I':'insertion_in_read','S':'soft_clip_reference','H':'soft_clip_read','=':'perfect_match','X':'mismatch'}
							read_seq_qual = sam_line[10]
							read_base_dict,CIGAR_operator_count_dict,barcode_base_string,barcode_match_string,barcode_qual_string,barcode_mismatches,backbone_base_string,backbone_match_string,backbone_qual_string,backbone_mismatches = SAM_parse_barcode(map_site,cigar,read_seq_qual,read_seq,map_qual,focal_region_start,focal_region_end,barcode_loc_dict,backbone_loc_dict)
							if len(read_base_dict) > 0 and barcode_base_string != '': #(read_nt,qual_score,mapping_qual,read_length,num,R_end_dist,nt_prop_location)
								barcode_region_reads_found += 1
								if barcode_mismatches == 0:
									barcode_reads_found += 1
									try:
										barcode_count_dict[barcode_base_string] += 1
									except:
										barcode_count_dict[barcode_base_string] = 1
								try:
									in_count = CIGAR_operator_count_dict['I']
								except:
									in_count = 0
								try:
									del_count = CIGAR_operator_count_dict['D']
								except:
									del_count = 0
								indel_count = in_count + del_count
								try:
									mismatch_count = CIGAR_operator_count_dict['X']
								except:
									mismatch_count = 0

								backbone_mismatch_count = len(backbone_match_string)-len(backbone_match_string.replace('X',''))
								barcode_outfile.write(barcode_base_string+"\t"+str(barcode_mismatches)+"\t"+str(backbone_mismatch_count)+"\t"+barcode_qual_string+"\t"+barcode_match_string+"\t"+backbone_base_string+"\t"+backbone_qual_string+"\t"+backbone_match_string+"\t"+str(in_count)+"\t"+str(del_count)+"\n")

			sam_infile.close()
			barcode_outfile.close()
			os.rename(temp_barcode_filename,barcode_filename_out)

			if dedupe_reads == True:
				os.remove(dedupe_samfilename)
			barcode_freq_array = []
			barcode_count_array = []
			for barcode_base_string in barcode_count_dict:
				count = barcode_count_dict[barcode_base_string]
				freq = float(count)/float(barcode_reads_found)
				barcode_count_array.append(count)
				barcode_freq_array.append(freq)
			if barcode_reads_found >=3:
				if len(barcode_count_dict) == 1:
					chao = 'na'
					shannon = ('na','na')
				else:
					chao = chao_1_richness(barcode_count_array)
					shannon = shannon_alpha(barcode_freq_array)
			else:
				chao = 'na'
				shannon = ('na','na')
			print(sampleID+" "+str(barcode_region_reads_found)+" "+str(barcode_reads_found)+" "+str(len(barcode_count_dict))+" "+str(shannon[0])+" "+str(shannon[1])+" "+str(chao))


####################### MAIN ##############################

bases = ['A','T','C','G']

continue_running = False
forward_read_suffix = ''
reverse_read_suffix = ''
suffix_info_file_path = temp_dir+"file_extensions.txt"
try:
	file_suffix_infile = open(suffix_info_file_path)
	for line in file_suffix_infile:
		line = line.strip().split("\t")
		if line[0] == "forward":
			forward_read_suffix = line[1]
		elif line[0] == "reverse":
			reverse_read_suffix = line[1]
		continue_running = True
except:
	forward_read_suffix = input('Read suffix info not found. Please enter the file suffix of the raw input reads to identify sampleIDs (or "exit" to stop):\n(e.g. for paired end reads that look like "GP01_day3_L001_R1.fastq", enter "_L001_R1.fastq")\n')
	if forward_read_suffix == "exit":
		sys.exit("Exiting.") 
	elif forward_read_suffix != '':
		continue_running = True

if continue_running == False:
	sys.exit("Exiting.")

## Create list of samples to process
file_list = [f for f in os.listdir(input_read_dir) if f.endswith(forward_read_suffix)]
sampleID_list = []
print_space_offset = 0
for files in file_list:
	sampleID = files.split(forward_read_suffix)[0]
	sampleID_list.append(sampleID)
	if len(sampleID) > print_space_offset:
		print_space_offset = len(sampleID)
sampleID_list = list(set(sampleID_list))

sample_ref_pair_list = []
for sampleID in sampleID_list:
	for ref_set in ref_set_list:
		primary_ref = ref_set.split("-")[0].split(".f")[0]
		if not os.path.exists(sub_dir+primary_ref):
			os.makedirs(sub_dir+primary_ref)
		if not os.path.exists(barcode_dir+primary_ref):
			os.makedirs(barcode_dir+primary_ref)
		if skip_qual == False:
			if not os.path.exists(qual_dir+primary_ref):
				os.makedirs(qual_dir+primary_ref)
		pair_tup = (sampleID,ref_set)
		sample_ref_pair_list.append(pair_tup)

processed_list = Parallel(n_jobs=num_cores)(delayed(main_pipeline)(sample_ref_tup,ref_set_list,input_read_dir,ref_dir,barcode_dir,trim_dir,map_dir,qual_dir,sam_dir,bam_dir,temp_dir,sub_dir,command_prefix,dedupe_reads,focal_region_start,focal_region_end,barcode_loc_dict,backbone_loc_dict) for sample_ref_tup in sample_ref_pair_list)

