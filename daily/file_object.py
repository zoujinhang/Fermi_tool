import os
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from Data_analysis import Clock
import Data_analysis.file as myfile
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


class Database(object):
	
	def __init__(self,databasepath,clock = None):
		'''
		
		:param databasepath:
		:param clock:
		'''
		self.topdir = databasepath
		self.detector = ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb','b0','b1']
		if clock is None:
			self.clock = Clock()
		else:
			self.clock = clock
		
	
	def input_pos_time(self,time_start,time_stop,format = None,scale='utc',astropyTime=False):
		'''
		
		:param time_start:
		:param time_stop:
		:param format:
		:param scale:
		'''
		if astropyTime==False:
			time_start = Time(time_start,format=format,scale=scale)
			time_stop =Time(time_stop,format=format,scale=scale)
		#mjd_in_h = np.arange(time_start.mjd,time_stop.mjd,1/24)
		#self.time_h_array = Time(mjd_in_h,format='mjd',scale='utc')
		mjd_in_day = np.arange(time_start.mjd,time_stop.mjd,1)
		self.time_day_array = Time(mjd_in_day,format='mjd',scale='utc')
		
		self.pos_table = self.get_poshist_data()

	def get_pos(self,t_array):
		'''
		
		:param t_array: met
		:return: pandas_table
		'''
		if self.pos_table is not None:
			t = np.array(t_array)

			c = {'T':t}
			pos_columns = self.pos_table.columns.values
			x = self.pos_table['SCLK_UTC']
			for y_name in pos_columns[1:]:
				y = self.pos_table[y_name]
				inter_f = interp1d(x,y,kind = 'cubic')
				c[y_name] = inter_f(t)
			return pd.DataFrame(c)
		else:
			print('The poshist file is missing and cannot provide posture information.')
			return None
		
	def get_detector_data(self,time_start,time_stop,format = None,scale='utc',astropyTime=False):
		'''
		
		:return: data {'n0':{'ch_E':pandas_table,            # columns CHANNEL,E_MIN,E_MAX
				'events':pandas_table},              #columns TIME, PHA
				...}
		'''
		if astropyTime==False:
			time_start = Time(time_start,format=format,scale=scale)
			time_stop =Time(time_stop,format=format,scale=scale)
		met_start = self.clock.utc_to_met(time_start,astropyTime = True)
		met_stop = self.clock.utc_to_met(time_stop,astropyTime = True)
		mjd_in_h = np.arange(time_start.mjd,time_stop.mjd,1/24)
		time_h_array = Time(mjd_in_h,format='mjd',scale='utc')
		date_time_arr = time_h_array.to_datetime()
		data = {}
		for deter in self.detector:#['n0','n1','n2','n3','n4','n5','n6','n7','n7','n8','n9','na','nb','b0','b1']
			data[deter] = {}
		
		for deter in self.detector:
			print('get',deter,'data.')
			#year,month,day
			data[deter]['ch_E'] = None
			data[deter]['events'] = None
			for date_t in date_time_arr:
				year = '%d' % date_t.year
				month = '%.2d' % date_t.month
				day = '%.2d' % date_t.day
				hour = '%.2d' % date_t.hour
				link = self.topdir + year + '/' + month + '/' + day + '/'
				name = 'glg_tte_'+deter+'_'+year[-2:]+month+day + '_'+hour+'z_v*'
				file = myfile.findfile(link,name)
				if len(file) > 0:
					file_name = file[0]
					if data[deter]['ch_E'] is None:
						hl = Table.read(link + file_name, hdu=1)
						data[deter]['ch_E'] = hl.to_pandas()
					if data[deter]['events'] is None:
						hl = Table.read(link + file_name, hdu=2)
						data[deter]['events'] = hl.to_pandas()
					else:
						hl = Table.read(link + file_name, hdu=2).to_pandas()
						data[deter]['events'].append(hl)
						#news = pd.concat([data[deter]['events'],hl])
						#news.drop_duplicates('TIME','first',inplace=True,ignore_index=True)
						#news = pd.merge(data[deter]['events'],hl,on = list(columns_)[0],how = 'outer')
						#data[deter]['events'] = news.sort_values(by = 'TIME')
				else:
					print('lost file:',name)
			data[deter]['events'].drop_duplicates('TIME','first',inplace=True,ignore_index=True)
			data[deter]['events'].sort_values(by = 'TIME',inplace=True,ignore_index=True)
			index_ = (data[deter]['events']['TIME']>=met_start)&(data[deter]['events']['TIME']<=met_stop)
			data[deter]['events'] = data[deter]['events'][index_]
	
		return data
				
	def get_poshist_data(self):
		'''
		
		:return:
		'''
		date_time_arr = self.time_day_array.to_datetime()
		timelist = []
		data = None
		for t_t in date_time_arr:
			year = '%d' % t_t.year
			month = '%.2d' % t_t.month
			day = '%.2d' % t_t.day
			link = self.topdir + year + '/' + month + '/' + day + '/'
			date = year+'-'+month+'-'+day
			timelist.append(date)
			name = 'glg_poshist_all_'+year[-2:]+month+day + '_v*'
			file = myfile.findfile(link,name)
			if len(file) > 0:
				file_name = file[0]
				if data is None:
					hl = Table.read(link + file_name,hdu = 1)
					data = hl.to_pandas()
				else:
					hl = Table.read(link + file_name,hdu = 1).to_pandas()
					data.append(hl)
					#columns_ = hl.columns.values
					#new = pd.concat([data,hl])
					#new.drop_duplicates('SCLK_UTC','first',inplace=True,ignore_index=True)
					#new = pd.merge(data,hl,on=columns_[:-1],how = 'outer')
					#new.sort_values(by = 'SCLK_UTC')
					#data = new
			else:
				print('The poshist file is missing.')
				print(name)
				return None
		data.drop_duplicates('SCLK_UTC','first',inplace=True,ignore_index=True)
		data.sort_values(by = 'SCLK_UTC')
		return data
	
	
	