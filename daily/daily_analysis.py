
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip,mad_std,bayesian_blocks
import operator
from Data_analysis.Baseline import TD_baseline
from Data_analysis.Bayesian_duration import background_correction,get_bayesian_duration
from scipy import stats
import pandas as pd
from astropy.coordinates import SkyCoord

class Sources(object):
	def __init__(self,positions=None,names=None,frame='icrs',unit = 'deg'):
		if positions is not None:
			if isinstance(positions,SkyCoord):
				self.positions = positions
			else:
				self.positions = SkyCoord(positions[0],positions[1],frame = frame,unit = unit)
			if names is not None:
				self.names = names
			else:
				try:
					self.names = np.arange(len(positions),dtype = int)
				except (TypeError):
					self.names = 0
	def input_with_deg(self,ra,dec,names = None,frame='icrs',unit = 'deg'):
		self.positions = SkyCoord(ra,dec,frame = frame,unit = unit)
		if names is not None:
			self.names = names
		else:
			try:
				self.names = np.arange(len(self.positions),dtype = int)
			except (TypeError):
				self.names = 0

def track(serch_result,geometry,sources):
	
	trig_a0 = serch_result['trig_0']
	#trig_a0.drop_duplicates('start','first',ignore_index=True)
	trig_a1 = serch_result['trig_1']
	#trig_a1.drop_duplicates('start','first',ignore_index=True)
	trig_all = serch_result['trig_all']
	sources_name = sources.names
	positions = sources.positions
	result ={'lc':serch_result['lc']}
	try:
		for i,pos in enumerate(positions):
			tl0 = track_one(trig_a0,geometry,pos)
			
			tl1 = track_one(trig_a1,geometry,pos)
			
			tl2 = track_one(trig_all,geometry,pos)
			
			result[sources_name[i]] = {'trig_0':trig_a0[tl0],
			                        'trig_1':trig_a1[tl1],
			                        'trig_all':trig_all[tl2]}
	except(TypeError):
		tl0 = track_one(trig_a0,geometry,positions)
		tl1 = track_one(trig_a1,geometry,positions)
		tl2 = track_one(trig_all,geometry,positions)
		result[sources_name] = {'trig_0':trig_a0[tl0],
			                'trig_1':trig_a1[tl1],
			                'trig_all':trig_all[tl2]}
				
	return 	result

def track_one(trig_a,geometry,pos):
	
	tool_ = []
	radius = geometry.radius
	t = trig_a['start'].values
	print(trig_a.shape)
	seq = geometry.get_separation_with_time(t, pos)
	try:
		ni_list = trig_a['overlap'].values
	except:
		ni_list = trig_a['detector'].values
		if seq is not None:
			for i,ni in enumerate(ni_list):
				seqi = seq.iloc[i]
				seq_values = seqi[ni]
				tool_.append(seq_values<=radius)
			
			return pd.Series(tool_,name='start',dtype=bool)
		else:
			print('seq is None!')
			return pd.Series([False]*trig_a.shape[0],name='start',dtype=bool)
	
	if seq is not None:
		for i,ni in enumerate(ni_list):
			n_ni = len(ni)
			seqi = seq.iloc[i]
			seq_values = seqi[ni].values
			ne = seq_values[seq_values<=radius]
			tool_.append(len(ne)>0.5*n_ni)
		return pd.Series(tool_,name='start',dtype=bool)
	else:
		print('seq is None!')
		return pd.Series([False]*trig_a.shape[0],name='start',dtype=bool)
		
		
		
def search_candidates(data,detectors,geometry):
	
	trig_data = {}
	lc = {}
	for deteri in detectors:
		ni = data[deteri]['events']
		t = ni['TIME'].values
		ni_c = analysis_one(t,binsize = 0.064,wt = 0.064,binsize_else=0.01,distinguish = 10,sigma = 3)
		trig_data[deteri] = ni_c
		lc[deteri] = {'lc':ni_c['lc'],'lc_bs':ni_c['lc_bs'],'sigma':ni_c['sigma']}
	c = trig_filrate(trig_data,geometry,detectors)
	c['lc'] = lc
	return c
	
def trig_filrate(trig_data,geometry,detectors):
	
	#detectors = geometry.detectors.name_list
	angle_overlap_ = geometry.detectors.get_angle_overlap()
	start_arr = np.array([])
	stop_arr = np.array([])
	wind_star = np.array([])
	wind_stop = np.array([])
	deter_name = []
	bayes = []
	for deteri in detectors:
		ni = trig_data[deteri]
		start_arr = np.concatenate((start_arr,ni['start'],ni['good_start']))
		stop_arr = np.concatenate((stop_arr,ni['stop'],ni['good_stop']))
		wind_star = np.concatenate((wind_star,ni['wind_start'],ni['good_wind_start']))
		wind_stop = np.concatenate((wind_stop,ni['wind_stop'],ni['good_wind_stop']))
		deter_name = deter_name + [deteri]*(len(ni['start'])+len(ni['good_start']))
		bayes = bayes + [0]*len(ni['start']) + [1]*len(ni['good_start'])
	
	trig_all = {'start':start_arr,'stop':stop_arr,'wind_start':wind_star,'wind_stop':wind_stop,'detector':deter_name,'bayes':bayes}
	trig_all = pd.DataFrame(trig_all)
	trig_all.sort_values(by='start',inplace = True,ignore_index=True)
	#print(tirg_all)
	index_0  = trig_all['bayes']==0
	tig_t = time_overlap(trig_all[index_0],n=3)
	tig_a0 = angle_overlap(tig_t,angle_overlap_,n = 3)
	index_1  = trig_all['bayes']==1
	tig_t = time_overlap(trig_all[index_1],n=2)
	tig_a1 = angle_overlap(tig_t,angle_overlap_,n = 2)
	c = {'trig_0':tig_a0,
	     'trig_1':tig_a1,
	     'trig_all':trig_all}
	return c
	
	
def angle_overlap(tig_all,angle_overlap,n = 3):
	'''
	
	:param tig_all:
	:param angle_overlap:
	:param n:
	:return:
	'''
	new_name_list = []
	new_start = []
	new_stop = []
	new_wind_start = []
	new_wind_stop = []
	new_overlap = []
	for i in range(tig_all.shape[0]):
		t_i = tig_all.iloc[i]
		ni_list = t_i['ni_list']
		v_ni_list = []
		for ni in ni_list:
			over_list = list(angle_overlap[ni])
			over_list.append(ni)
			v_ni_list = v_ni_list+over_list
		ni_n = list(set(v_ni_list))     #Remove duplicates
		ni_n = np.array(ni_n)
		ni_nm =[]
		v_ni_list = np.array(v_ni_list)
		for ni in ni_n:
			index_n = np.where(v_ni_list == ni)[0]
			nn = v_ni_list[index_n].size
			ni_nm.append(nn)
		ni_nm = np.array(ni_nm)
		index_ = np.where(ni_nm>1)[0]
		ni_n = ni_n[index_]                         #get overlap
		ni_over_set = set(ni_n)
		ni_set = set(ni_list)
		ni_union_set = list(ni_set & ni_over_set)   #get overlap
		num = len(ni_union_set)
		if num>n:
			new_name_list = new_name_list+ni_union_set
			new_start = new_start + [t_i['start']]*num
			new_stop = new_stop + [t_i['stop']]*num
			new_wind_start = new_wind_start + [t_i['wind_start']]*num
			new_wind_stop = new_wind_stop + [t_i['wind_stop']]*num
			new_overlap = new_overlap + [ni_union_set]*num
	c = {'start':np.array(new_start),
	     'stop':np.array(new_stop),
	     'wind_start':np.array(new_wind_start),
	     'wind_stop':np.array(new_wind_stop),
	     'detector':new_name_list,
	     'overlap':new_overlap}
	return pd.DataFrame(c)
		
def time_overlap(tig_all,n = 3):
	'''
	
	:param tig_all:
	:param n:
	:return:
	'''
	new_start = []
	new_stop = []
	new_name = []
	new_wind_start = []
	new_wind_stop = []
	ni_list = []
	t_over_list = []
	for i in range(tig_all.shape[0]):
		t_i = tig_all.iloc[i]
		ni = t_i['detector']
		start = t_i['start']
		stop = t_i['stop']
	
		if i == 0 :
			t_over_list.append([start,stop])
			ni_list.append(ni)
		else:
			trun_n = 0
			for start0,stop0 in t_over_list:
				if ((start<stop0)and(stop>start0)):#time overlap.
					trun_n = trun_n+1
			if len(t_over_list) == trun_n:
				t_over_list.append([start,stop])
				ni_list.append(ni)
			elif trun_n>=2:
				t_over_list.append([start,stop])
				ni_list.append(ni)
			elif trun_n == 0:
				if len(t_over_list) >= n:
					t_over_array = np.array(t_over_list)
					t_max = t_over_array.max()
					t_min = t_over_array.min()
					t_during = t_max - t_min
					wind_during_harf = 0.5*t_during/0.62
					if wind_during_harf<2.5:
						wind_during_harf=2.5
					t_center = 0.5*(t_min+t_max)
					new_start.append(t_min)
					new_stop.append(t_max)
					new_wind_start.append(t_center-wind_during_harf)
					new_wind_stop.append(t_center+wind_during_harf)
					new_name.append(list(set(ni_list)))
					t_over_list = [[start,stop]]
					ni_list = [ni]
				else:
					t_over_list = [[start,stop]]
					ni_list = [ni]
			else:
				if len(t_over_list) >= n:
					t_over_array = np.array(t_over_list)
					t_max = t_over_array.max()
					t_min = t_over_array.min()
					t_during = t_max - t_min
					wind_during_harf = 0.5 * t_during / 0.62
					if wind_during_harf < 2.5:
						wind_during_harf = 2.5
					t_center = 0.5 * (t_min + t_max)
					new_start.append(t_min)
					new_stop.append(t_max)
					new_wind_start.append(t_center - wind_during_harf)
					new_wind_stop.append(t_center + wind_during_harf)
					new_name.append(list(set(ni_list)))
					t_over_list = [[start,stop]]
					ni_list = [ni]
				else:
					t_over_list = [[start,stop]]
					ni_list = [ni]
	c = {'start':np.array(new_start),
	     'stop':np.array(new_stop),
	     'wind_start':np.array(new_wind_start),
	     'wind_stop':np.array(new_wind_stop),
	     'ni_list':new_name}
	return pd.DataFrame(c)
	

def analysis_one(t,binsize = 0.064,wt = 0.064,binsize_else = 0.01,distinguish=3,sigma = 3):
	'''
	
	:param t:
	:param binsize:
	:param wt:
	:param binsize_else:
	:param distinguish:
	:param sigma:
	:return:
	'''
	lc_list,t_index_list = get_light_curve_list(t,binsize = binsize,wt = wt)
	good_wind_start = []
	good_wind_stop = []
	good_start = []
	good_stop = []
	wind_start = []
	wind_stop = []
	start = []
	stop = []
	lc_bs_list = []
	sigma_list = []
	for lc in lc_list:
		lc_t,lc_rate = lc
		lc_t,lc_cs,lc_bs = TD_baseline(lc_t,lc_rate)
		lc_bs_list.append(lc_bs)
		mask = sigma_clip(lc_cs,sigma=5,maxiters=5,stdfunc=mad_std).mask
		myfilter = list(map(operator.not_, mask))
		lc_median_part = lc_cs[myfilter]
		loc,scale = stats.norm.fit(lc_median_part)
		sigma_list.append(scale)
		index_ = np.where(lc_cs>loc+sigma*scale)[0]
		if len(index_)>0:
			lc_t_list = []
			lc_cs_new = lc_cs + lc_bs.mean()
			i_list = get_subsection_index(index_,binsize,distinguish)
			for i in i_list:
				lc_ti = lc_t[i]
				if len(lc_ti) <5:
					lc_t_list.append([lc_ti.min()-2*binsize,lc_ti.max()+2*binsize])
				else:
					lc_t_list.append([lc_ti.min(),lc_ti.max()])
			#print('lc_t_list',lc_t_list)
			for ind,lc_ti in enumerate(lc_t_list):
				lc_dt = lc_ti[-1]-lc_ti[0]
				#print('lc_dt',lc_dt)
				if lc_dt<=2:
					add_t = 15
				elif lc_dt<=50:
					add_t = 80
				else:
					add_t = 100
				range_t_min = lc_ti[0]-add_t
				if ind == 0:
					if range_t_min < lc_t.min():
						range_t_min = lc_t.min()
				else:
					if range_t_min < lc_t_list[ind-1][-1]:
						range_t_min = lc_t_list[ind-1][-1]
				range_t_max = lc_ti[-1] + add_t
				if ind == len(lc_t_list)-1:
					if range_t_max >lc_t.max():
						range_t_max = lc_t.max()
				else:
					if range_t_max >lc_t_list[ind+1][0]:
						range_t_max = lc_t_list[ind+1][0]
						
				#print('range_t_min',range_t_min)
				#print('range_t_max',range_t_max)
				#print('d_range_t',range_t_max-range_t_min)
				t_index = np.where((lc_t>=range_t_min)&(lc_t<=range_t_max))[0]
				lt_t = lc_t[t_index]
				#print('lt_t',lt_t)
				#lt_rate = lc_rate[t_index]
				lt_rate_new = lc_cs_new[t_index]
				edges = bayesian_blocks(lt_t,np.round(lt_rate_new*binsize),fitness='events',p0 = 0.05,gamma = np.exp(-5))
				if len(edges)>=4:
					result = background_correction(lt_t,lt_rate_new,edges,degree = 6)
					startedges,stopedges = get_bayesian_duration(result,sigma = 3)
					if startedges.size == stopedges.size:
						if startedges.size >0:
							good_wind_start = good_wind_start+[range_t_min]*startedges.size
							good_wind_stop = good_wind_stop+[range_t_max]*startedges.size
							good_start = good_start + list(startedges)
							good_stop = good_stop + list(stopedges)
						else:
							wind_start.append(range_t_min)
							wind_stop.append(range_t_max)
							start.append(lc_ti[0])
							stop.append(lc_ti[-1])
					else:
						wind_start.append(range_t_min)
						wind_stop.append(range_t_max)
						start.append(lc_ti[0])
						stop.append(lc_ti[-1])
				else:
					if add_t == 10:#Temporary abandonment
						print('work in events!')
						t_index = np.where((t >= range_t_min) & (t <= range_t_max))[0]
						t_in = t[t_index]
						edges = bayesian_blocks(t_in, fitness='events', p0=0.01)
						gg = np.around(edges / binsize_else)
						gg = np.unique(gg)
						gg = np.sort(gg)
						edges = gg * binsize_else
						if len(edges) >= 4:
							
							lt_bins = np.arange(t_in.min(), t_in.max(),
							                    binsize_else)
							lt_bin_n = np.histogram(t_in, bins=lt_bins)[0]
							lt_bin_c = 0.5 * (lt_bins[1:] + lt_bins[:-1])
							lt_t, lt_cs, lt_bs = TD_baseline(lt_bin_c,
							                                 lt_bin_n / binsize_else)
							lt_rate = lt_cs + lt_bs.mean()
							result = background_correction(lt_t, lt_rate, edges,
							                               degree=6)
							startedges, stopedges = get_bayesian_duration(result,
							                                              sigma=3)
							if startedges.size == stopedges.size:
								if startedges.size > 0:
									good_wind_start = good_wind_start + [
										range_t_min] * startedges.size
									good_wind_stop = good_wind_stop + [
										range_t_max] * startedges.size
									good_start = good_start + list(startedges)
									good_stop = good_stop + list(stopedges)
								else:
									wind_start.append(range_t_min)
									wind_stop.append(range_t_max)
									start.append(lc_ti[0])
									stop.append(lc_ti[-1])
							else:
								wind_start.append(range_t_min)
								wind_stop.append(range_t_max)
								start.append(lc_ti[0])
								stop.append(lc_ti[-1])
						else:
							wind_start.append(range_t_min)
							wind_stop.append(range_t_max)
							start.append(lc_ti[0])
							stop.append(lc_ti[-1])
					else:
						wind_start.append(range_t_min)
						wind_stop.append(range_t_max)
						start.append(lc_ti[0])
						stop.append(lc_ti[-1])
	c = {
		'lc':lc_list,
		'lc_bs':lc_bs_list,
		'sigma':sigma_list,
		'good_wind_start':np.array(good_wind_start),
		'good_wind_stop':np.array(good_wind_stop),
		'good_start':np.array(good_start),
		'good_stop':np.array(good_stop),
		'wind_start':np.array(wind_start),
		'wind_stop':np.array(wind_stop),
		'start':np.array(start),
		'stop':np.array(stop)
	}
	return c
	
	
		
def get_subsection_index(index,binsize,distinguish=3):
	'''
	
	:param index:
	:param binsize:
	:param distinguish:
	:return:
	'''
	if len(index) == 1:
		return [index]
	else:
		index_list = []
		index_one = []
		d_index = index[1:]-index[:-1]
		d_index = np.concatenate(([0],d_index))
		for id,value in enumerate(d_index):
			if value*binsize <= distinguish:
				index_one.append(index[id])
			else:
				index_list.append(np.array(index_one))
				index_one = [index[id]]
		index_list.append(np.array(index_one))
		return index_list
	

def get_light_curve_list(t,binsize = 0.064,wt = 0.064):
	t_index_list = data_chack(t,wt = wt)
	t_index_list_new = []
	lc_list = []
	for t_index in t_index_list:
		t_i = t[t_index]
		min_ = t_i.min()
		max_ = t_i.max()
		if (max_ - min_) > 10*binsize:
			bins = np.arange(min_,max_,binsize)
			bin_n = np.histogram(t_i,bins=bins)[0]
			bin_c = 0.5*(bins[1:]+bins[:-1])
			lc_list.append([bin_c,bin_n/binsize])
			t_index_list_new.append(t_index)
	#print(t_index_list)
	return lc_list,t_index_list_new
	

def data_chack(t,wt = 0.064):
	#print('len t',len(t))
	
	dt = t[1:]-t[:-1]
	dt = np.concatenate((dt[:1],dt))
	#index_d = np.where(dt>=wt)[0]
	#if len(index_d)>0
	#print('len dt',len(dt))
	#print('dt',dt)
	index_list = []
	index_ = []
	not_first = False
	#print('all:',dt.size)
	for index,dti in enumerate(dt):
		
		if dti < wt:
			not_first = True
			index_.append(index)
		else:
			#print('chack!')
			if not_first:
				#print('add',len(index_))
				index_list.append(np.array(index_))
				index_ = [index]
			
	#print('add',len(index_))
	index_list.append(np.array(index_))
	return index_list


























