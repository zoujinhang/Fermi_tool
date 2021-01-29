#import os
from astropy.io import fits
#from astropy.table import Table
from astropy.time import Time
from Data_analysis import Clock
import Data_analysis.file as myfile
import numpy as np
import pandas as pd
#from scipy.interpolate import interp1d

def time_slic(time_start,time_stop,format = None,scale='utc',H = 3):
	
	if isinstance(time_start,Time) ==False:
		time_start = Time(time_start, format=format, scale=scale)
	if isinstance(time_stop,Time) ==False:
		time_stop =Time(time_stop,format=format,scale=scale)
	re = pd.date_range(time_start.fits,time_stop.fits,freq = str(H)+'H')
	re = pd.DataFrame({'time':re})
	stp = pd.date_range(time_stop.fits,freq = 'H',periods=1)
	stp = pd.DataFrame({'time':stp})
	
	t_ar = pd.concat([re,stp])
	t_ar.drop_duplicates('time','first',inplace=True,ignore_index=True)
	t_ar.sort_values(by = 'time',inplace=True,ignore_index=True)
	time_list = t_ar['time'].values
	t_ar1 = pd.DataFrame({'time_start':time_list[:-1],'time_stop':time_list[1:]})
	return t_ar1.values
	
	
class Database(object):

	def __init__(self,databasepath,detectors = None,clock = None):

		self.topdir = databasepath
		if detectors is not None:
			self.detector = detectors
		else:
			self.detector = ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb','b0','b1']

		if clock is None:

			self.clock = Clock()

		else:
			self.clock = clock

	def get_detector_data(self,time_start,time_stop,format = None,scale='utc'):

		'''

		:return: data {'n0':{'ch_E':pandas_table,            # columns CHANNEL,E_MIN,E_MAX
				'events':pandas_table},              #columns TIME, PHA
				...}
		'''
		if isinstance(time_start,Time) ==False:
			time_start = Time(time_start, format=format, scale=scale)
		if isinstance(time_stop,Time) ==False:
			time_stop =Time(time_stop,format=format,scale=scale)
		met_start = self.clock.utc_to_met(time_start)-1
		met_stop = self.clock.utc_to_met(time_stop)+1
		t_start = pd.Timestamp(time_start.fits)
		year = '%d' % t_start.year
		month = '%.2d' % t_start.month
		day = '%.2d' % t_start.day
		hour = '%.2d' % t_start.hour

		new_start  = year + '-' + month + '-' + day + 'T' + hour + ':00:00'

		t_stop = pd.Timestamp(time_stop.fits)

		year = '%d' % t_stop.year
		month = '%.2d' % t_stop.month
		day = '%.2d' % t_stop.day
		hour = '%.2d' % t_stop.hour

		new_stop = year + '-' + month + '-' + day + 'T' + hour + ':00:00'

		date_time_arr = pd.date_range(new_start,new_stop,freq = 'H')
		print(date_time_arr)
		data_c = {}
		for deter in self.detector:#['n0','n1','n2','n3','n4','n5','n6','n7','n7','n8','n9','na','nb','b0','b1']
			data_c[deter] = {}

		for deter in self.detector:
			data_c[deter]['ch_E'] = None
			data_c[deter]['events'] = None
			t = []
			ch = []
			for date_i in range(date_time_arr.shape[0]):
				date_t = date_time_arr[date_i]
				year = '%d' % date_t.year
				month = '%.2d' % date_t.month
				day = '%.2d' % date_t.day
				hour = '%.2d' % date_t.hour
				utc_start_time = year + '-' + month + '-' + day + 'T' + hour + ':00:00'
				met_starti = self.clock.utc_to_met(utc_start_time)
				met_stopi = met_starti + 3600.0
				link = self.topdir + year + '/' + month + '/' + day + '/'
				name = 'glg_tte_'+deter+'_'+year[-2:]+month+day + '_'+hour+'z_v*'
				file = myfile.findfile(link,name)
				if len(file) > 0:
					file_name = file[0]
					hl = fits.open(link + file_name)
					if data_c[deter]['ch_E'] is None:

						data = hl[1].data
						data_c[deter]['ch_E'] = pd.DataFrame({'CHANNEL':np.array(data['CHANNEL'],dtype=np.int16),
										    'E_MIN':np.array(data['E_MIN'],dtype = np.float),
										    'E_MAX':np.array(data['E_MAX'],dtype = np.float)})
					data = hl[2].data
					indexi = np.where((data['TIME']>=met_starti)&(data['TIME']<=met_stopi))[0]
					t.append(np.array(data['TIME'][indexi],dtype = np.float))
					ch.append(np.array(data['PHA'][indexi],dtype = np.int16))
					hl.close()
				else:
					print('lost file:',name)
			if data_c[deter]['ch_E'] is not None:
				t = np.concatenate(t)
				ch = np.concatenate(ch)
				sort_index = np.argsort(t)
				t = t[sort_index]
				ch = ch[sort_index]
				_,un_index = np.unique(t,return_index=True)
				t = t[un_index]
				ch = ch[un_index]
				t_index = np.where((t>=met_start)&(t<=met_stop))[0]
				data_c[deter]['events'] = pd.DataFrame({'TIME':t[t_index],
								      'PHA':ch[t_index]})
			else:
				data_c[deter]['events'] = None
		return data_c


	def get_poshist_data(self,time_start,time_stop,format = None,scale='utc'):

		if isinstance(time_start,Time) == False:

			time_start = Time(time_start, format=format, scale=scale)

		if isinstance(time_stop,Time) == False:

			time_stop =Time(time_stop,format=format,scale=scale)

		met_start = self.clock.utc_to_met(time_start)-1
		met_stop = self.clock.utc_to_met(time_stop)+1
		t_start = pd.Timestamp(time_start.fits)
		year = '%d' % t_start.year
		month = '%.2d' % t_start.month
		day = '%.2d' % t_start.day
		#hour = '%.2d' % t_start.hour

		new_start = year + '-' + month + '-' + day + 'T00:00:00'

		t_stop = pd.Timestamp(time_stop.fits)

		year = '%d' % t_stop.year
		month = '%.2d' % t_stop.month
		day = '%.2d' % t_stop.day
		#hour = '%.2d' % t_stop.hour

		new_stop = year + '-' + month + '-' + day + 'T00:00:00'
		date_time_arr = pd.date_range(new_start,new_stop,freq = 'D')
		t = []
		q1,q2,q3,q4 = [],[],[],[]
		x,y,z = [],[],[]
		lat,lon = [],[]
		col_name = ['SCLK_UTC','QSJ_1','QSJ_2','QSJ_3','QSJ_4','POS_X','POS_Y','POS_Z','SC_LAT','SC_LON']
		for t_i in range(date_time_arr.shape[0]):

			year = '%d' % date_time_arr[t_i].year
			month = '%.2d' % date_time_arr[t_i].month
			day = '%.2d' % date_time_arr[t_i].day
			link = self.topdir + year + '/' + month + '/' + day + '/'
			name = 'glg_poshist_all_'+year[-2:]+month+day + '_v*'
			file = myfile.findfile(link,name)
			if len(file) > 0:
				file_name = file[0]
				hl = fits.open(link + file_name)
				data = hl[1].data
				t.append(data['SCLK_UTC'])
				q1.append(data['QSJ_1'])
				q2.append(data['QSJ_2'])
				q3.append(data['QSJ_3'])
				q4.append(data['QSJ_4'])
				x.append(data['POS_X'])
				y.append(data['POS_Y'])
				z.append(data['POS_Z'])
				lat.append(data['SC_LAT'])
				lon.append(data['SC_LON'])
				hl.close()
			else:
				print('warting! the file under line is lost!')
				print(name)
		if len(t)>0:
			t = np.concatenate(t)
			q1,q2,q3,q4 = np.concatenate(q1),np.concatenate(q2),np.concatenate(q3),np.concatenate(q4)
			x,y,z = np.concatenate(x),np.concatenate(y),np.concatenate(z)
			lat,lon = np.concatenate(lat),np.concatenate(lon)
			inde_sort = np.argsort(t)
			new_data = np.vstack([t,q1,q2,q3,q4,x,y,z,lat,lon]).T
			new_data = new_data[inde_sort]
			_,uni_index = np.unique(new_data[:,0],return_index=True)
			new_data = new_data[uni_index]
			index_ = np.where((new_data[:,0]>=met_start-2)&(new_data[:,0]<=met_stop+2))[0]
			new_data = new_data[index_]
			d_data = pd.DataFrame(data=new_data,dtype=np.float,columns=col_name)
			return d_data
		else:
			return None

