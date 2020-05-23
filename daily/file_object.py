import os
from astropy.io import fits
from astropy.time import Time
from Data_analysis import Clock
import Data_analysis.file as myfile
import numpy as np

class Database(object):
	
	def __init__(self,databasepath,clock = None):
		'''
		
		:param databasepath:
		:param clock:
		'''
		self.topdir = databasepath
		self.detector = ['n0','n1','n2','n3','n4','n5','n6','n7','n7','n8','n9','na','nb','b0','b1']
		if clock is None:
			self.clock = Clock()
		else:
			self.clock = clock
	
	def get_files(self,time_start,time_stop,format = None,scale='utc'):
		'''
		
		:param time_start:
		:param time_stop:
		:param format:
		:param scale:
		:return:
		'''
		time_start = Time(time_start,format=format,scale=scale)
		time_stop =Time(time_stop,format=format,scale=scale)
		mjd_in_h = np.arange(time_start.mjd,time_stop.mjd+1/24,1/24)
		time_h_array = Time(mjd_in_h,format='mjd',scale='utc')
		mjd_in_day = np.arange(time_start.mjd,time_stop.mjd+1,1)
		time_day_array = Time(mjd_in_day,format='mjd',scale='utc')
		t_list_h,data_dete = self.get_detector_file(time_h_array)
		t_list_day,data_pos = self.get_poshist_file(time_day_array)
		retr = {
			'pos':[t_list_day,data_pos],
		        'det':[t_list_h,data_dete]
		        }
		return retr
		
		
	def get_detector_file(self,time_array):
		'''
		
		:param time_array:
		:return:
		'''
		
		date_time_arr = time_array.to_datetime()
		timelist = []
		data = {}
		for deter in self.detector:
			data[deter] = {}
		for date_t in date_time_arr:
			#year,month,day
			year = '%d' % date_t.year
			month = '%2d' % date_t.month
			day = '%2d' % date_t.day
			hour = '%2d' % date_t.hour
			link = self.topdir + year + '/' + month + '/' + day + '/'
			date = year+'-'+month+'-'+day+'T'+hour
			timelist.append(date)
			for deter in self.detector:
				name = 'glg_tte_'+deter+'_'+year[-2:]+month+day + '_'+hour+'z_v*'
				file = myfile.findfile(link,name)
				if len(file) > 0:
					file_name = file[0]
					hl = fits.open(link + file_name)
					data[deter][date] = hl
				else:
					data[deter][date] = None
		return timelist,data
				
	def get_poshist_file(self,time_array):
		'''
		
		:param time_array:
		:return:
		'''
		date_time_arr = time_array.to_datetime()
		timelist = []
		data = {}
		for t_t in date_time_arr:
			year = '%d' % t_t.year
			month = '%2d' % t_t.month
			day = '%2d' % t_t.day
			link = self.topdir + year + '/' + month + '/' + day + '/'
			date = year+'-'+month+'-'+day
			timelist.append(date)
			name = 'glg_poshish_all_'+year[-2:]+month+day + '_v*'
			file = myfile.findfile(link,name)
			if len(file) > 0:
				file_name = file[0]
				hl = fits.open(link + file_name)
				data[date] = hl
			else:
				data[date] = None
		return timelist,data
	
	
	