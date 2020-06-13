
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

class Plot_track(object):
	
	def __init__(self,result,detector,geometry):
		'''
		
		:param result:
		:param detector: detector name list
		:param geometry:
		'''
		self.result = result
		self.detector = detector
		self.geometry = geometry
		self.clock = geometry.Time_transition
		self.t_r_cmap = ListedColormap(['k', '#d3d7d4'])
		self.t_r_norm = BoundaryNorm([0, self.geometry.radius,181], self.t_r_cmap.N)
		self.t_b_cmap = ListedColormap(['#f47920', '#ffce7b'])
		self.t_b_norm = BoundaryNorm([0, self.geometry.radius,181], self.t_b_cmap.N)
		
		
	def plot_one_source(self,sn,positions,axs):
		'''
		
		:param sn:
		:param positions:
		:param axs:
		:return:
		'''
		
		
		tirg_data = self.result['lc']
		sn_serch_result = self.result[sn]
		sn_trig_0 = sn_serch_result['trig_0']
		sn_trig_1 = sn_serch_result['trig_1']
		sn_trig_all = sn_serch_result['trig_all']
		for ind,deteri in enumerate(self.detector):
			n0_c = tirg_data[deteri]
			det_trig_0 = sn_trig_0[sn_trig_0['detector'] == deteri]
			det_trig_1 = sn_trig_1[sn_trig_1['detector'] == deteri]
			det_trig_all = sn_trig_all[sn_trig_all['detector'] == deteri]
			lc_list = n0_c['lc']
			lc_bs_list = n0_c['lc_bs']
			t_min = lc_list[0][0][0]
			utc_start = self.clock.met_to_utc(t_min).fits
			t_max = lc_list[-1][0][-1]
			start = det_trig_all['start'].values
			stop  = det_trig_all['stop'].values
			bayes = det_trig_all['bayes'].values
			for i in range(len(start)):
				if bayes[i] == 0:
					axs[ind * 2].axvline(x=start[i] -t_min, color='#4f5555', alpha=0.2)
					axs[ind * 2].axvline(x=stop[i] - t_min, color='#4f5555', alpha=0.2)
				else:
					axs[ind * 2].axvline(x=start[i] - t_min, color='#c37e00', alpha=0.2)
					axs[ind * 2].axvline(x=stop[i] - t_min, color='#c37e00', alpha=0.2)
			
			start = det_trig_0['start'].values
			stop = det_trig_0['stop'].values
			for i in range(len(start)):
				axs[ind * 2].axvline(x=start[i] - t_min, color='#ea66a6')
				axs[ind * 2].axvline(x=stop[i] - t_min, color='#009ad6')
			
			start = det_trig_1['start'].values
			stop = det_trig_1['stop'].values
			for i in range(len(start)):
				axs[ind * 2].axvline(x=start[i] - t_min, color='#aa2116')
				axs[ind * 2].axvline(x=stop[i] - t_min, color='#225a1f')

			
			axs[ind*2].plot(0,0,color = 'k',label = deteri+' lightcurve')
			axs[ind*2].plot(0,0,color = '#f47920',label = 'background')
			lc_ratemax = 0
			for indx,lc in enumerate(lc_list):
				lc_t,lc_rate = lc
				lcr_max = lc_rate.max()
				if lcr_max>lc_ratemax:
					lc_ratemax = lcr_max
					
				lc_bs = lc_bs_list[indx]
				
				lc_tc = 0.5*(lc_t[:-1]+lc_t[1:])
				
				t_r_points = np.array([lc_t-t_min, lc_rate]).T.reshape(-1, 1, 2)
				t_r_segments = np.concatenate([t_r_points[:-1], t_r_points[1:]], axis=1)
				
				t_b_points = np.array([lc_t-t_min, lc_bs]).T.reshape(-1, 1, 2)
				
				t_b_segments = np.concatenate([t_b_points[:-1], t_b_points[1:]], axis=1)
				
				seq = self.geometry.get_separation_with_time(lc_tc,positions)
				seq = seq[deteri].values
				
				t_r = LineCollection(t_r_segments, cmap=self.t_r_cmap, norm=self.t_r_norm)
				t_r.set_array(seq)
				t_r.set_linewidth(2)
				#lc.set_linewidth(2)
				axs[ind*2].add_collection(t_r)

				t_b = LineCollection(t_b_segments, cmap=self.t_b_cmap, norm=self.t_b_norm)
				t_b.set_array(seq)
				t_b.set_linewidth(2)
				axs[ind*2].add_collection(t_b)
				axs[ind*2+1].plot(lc_tc-t_min,seq,color = 'k')
			axs[ind*2].set_xticks([])
			axs[ind*2].set_xlim(0,t_max-t_min)
			axs[ind*2].set_ylim(0,lc_ratemax/0.8)
			axs[ind*2].set_ylabel('Count rate')
			axs[ind*2].legend(loc = 'upper left')
			axs[ind*2+1].set_xlim(0,t_max-t_min)
			axs[ind*2+1].set_ylim(0,180)
			axs[ind*2+1].set_ylabel('Separation degree')
			axs[ind*2+1].set_xlabel('Start at ' + str(utc_start)+' (s)')
			#axs[ind*2+1].legend()
	
	def plot_one_detector(self,sn,detector,positions,axs):
		'''
		
		:param sn:
		:param detector:
		:param positions:
		:param axs:
		:return:
		'''
		tirg_data = self.result['lc']
		sn_serch_result = self.result[sn]
		sn_trig_0 = sn_serch_result['trig_0']
		sn_trig_1 = sn_serch_result['trig_1']
		sn_trig_all = sn_serch_result['trig_all']
		n0_c = tirg_data[detector]
		lc_list = n0_c['lc']
		lc_bs_list = n0_c['lc_bs']
		t_min = lc_list[0][0][0]
		utc_start = self.clock.met_to_utc(t_min).fits
		t_max = lc_list[-1][0][-1]
		
		det_trig_0 = sn_trig_0[sn_trig_0['detector'] == detector]
		det_trig_1 = sn_trig_1[sn_trig_1['detector'] == detector]
		det_trig_all = sn_trig_all[sn_trig_all['detector'] == detector]
		
		start = det_trig_all['start'].values
		stop  = det_trig_all['stop'].values
		bayes = det_trig_all['bayes'].values
		
		for i in range(len(start)):
			if bayes[i] == 0:
				axs[0].axvline(x=start[i] - t_min, color='#4f5555', alpha=0.2)
				axs[0].axvline(x=stop[i] - t_min, color='#4f5555', alpha=0.2)
			else:
				axs[0].axvline(x=start[i] - t_min, color='#c37e00', alpha=0.2)
				axs[0].axvline(x=stop[i] - t_min, color='#c37e00', alpha=0.2)
			
		start = det_trig_0['start'].values
		stop = det_trig_0['stop'].values
		for i in range(len(start)):
			axs[0].axvline(x=start[i] - t_min, color='#ea66a6')
			axs[0].axvline(x=stop[i] - t_min, color='#009ad6')
		
		start = det_trig_1['start'].values
		stop = det_trig_1['stop'].values
		for i in range(len(start)):
			axs[0].axvline(x=start[i] - t_min, color='#aa2116')
			axs[0].axvline(x=stop[i] - t_min, color='#225a1f')

		
		axs[0].plot(0,0,color = 'k',label = detector+' lightcurve')
		axs[0].plot(0,0,color = '#f47920',label = 'background')
		lc_ratemax = 0
		for indx,lc in enumerate(lc_list):
			lc_t,lc_rate = lc
			lcr_max = lc_rate.max()
			if lcr_max>lc_ratemax:
				lc_ratemax = lcr_max
			lc_bs = lc_bs_list[indx]
			lc_tc = 0.5 * (lc_t[:-1] + lc_t[1:])
			
			t_r_points = np.array([lc_t - t_min, lc_rate]).T.reshape(-1, 1, 2)
			t_r_segments = np.concatenate([t_r_points[:-1], t_r_points[1:]], axis=1)
			
			t_b_points = np.array([lc_t - t_min, lc_bs]).T.reshape(-1, 1, 2)
			t_b_segments = np.concatenate([t_b_points[:-1], t_b_points[1:]], axis=1)
			
			seq = self.geometry.get_separation_with_time(lc_tc, positions)
			seq = seq[detector].values
			
			t_r = LineCollection(t_r_segments, cmap=self.t_r_cmap, norm=self.t_r_norm)
			t_r.set_array(seq)
			t_r.set_linewidth(2)
			axs[0].add_collection(t_r)
			
			t_b = LineCollection(t_b_segments, cmap=self.t_b_cmap, norm=self.t_b_norm)
			t_b.set_array(seq)
			t_b.set_linewidth(2)
			axs[0].add_collection(t_b)
			
			axs[1].plot(lc_tc - t_min, seq, color='k')
		axs[0].set_xticks([])
		axs[0].set_xlim(0,t_max-t_min)
		axs[0].set_ylim(0,lc_ratemax/0.8)
		axs[0].set_ylabel('Count rate')
		axs[0].legend(loc = 'upper left')
		axs[1].set_xlim(0,t_max-t_min)
		axs[1].set_ylim(0,180)
		axs[1].set_ylabel('Separation degree')
		axs[1].set_xlabel('Start at ' + str(utc_start)+' (s)')



