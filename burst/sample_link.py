import os
import numpy as np
import Data_analysis.file as myfile
from astropy.io import fits


def get_sample_dir_list(yearlist, databaselink):
	
	'''
	
	:param yearlist:
	:param databaselink:
	:return:
	'''
	sample_dir_list = []
	for year in yearlist:
		topdir = databaselink + str(year) + '/'
		dirlist1 = os.listdir(topdir)
		dirlist1 = np.sort(dirlist1)
		for dirl in dirlist1:
			if os.path.isdir(topdir + dirl):
				sample_dir_list.append([topdir + dirl + '/', dirl])
	return sample_dir_list

def get_file(input_list, NaI, BGO):
	
	'''
	
	:param input_list:
	:param NaI:
	:param BGO:
	:return:
	'''
	
	sample = input_list[1]
	sampledir = input_list[0]
	data = {}
	loc_name_list = myfile.findfile(sampledir, 'glg_locprob_all_' + sample + '_v*')
	if len(loc_name_list) >= 1:
		data['loc'] = fits.open(sampledir + loc_name_list[0])
	else:
		loc_name_list = myfile.findfile(sampledir, 'glg_tcat_all_' + sample + '_v*')
		if len(loc_name_list) >= 1:
			data['loc'] = fits.open(sampledir + loc_name_list[0])
		else:
			data['loc'] = None
	trigdat_name_list = myfile.findfile(sampledir, 'glg_trigdat_all_' + sample + '_v*')
	if len(trigdat_name_list) >= 1:
		data['trigdat'] = fits.open(sampledir + trigdat_name_list[0])
	else:
		data['trigdat'] = None
	for ni in NaI:
		ni_list = myfile.findfile(sampledir, 'glg_tte_' + ni + '_' + sample + '_v*')
		if len(ni_list) >= 1:
			# print(ni)
			data[ni] = fits.open(sampledir + ni_list[0])
		else:
			data[ni] = None
	for bi in BGO:
		bi_list = myfile.findfile(sampledir, 'glg_tte_' + bi + '_' + sample + '_v*')
		if len(bi_list) >= 1:
			# print(bi)
			data[bi] = fits.open(sampledir + bi_list[0])
		else:
			data[bi] = None
	return data

