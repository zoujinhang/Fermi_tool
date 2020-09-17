
import numpy as np
from Data_analysis import TD_baseline
from Data_analysis import overlap_curve


def get_spectrum(t, ch, ch_n, bin_arr = None, edges=None, bg_dt=1.0):
	data = []
	data_err = []
	bins = np.arange(t[0], t[-1], bg_dt)
	for chi in ch_n:
		index = np.where(ch == chi)
		t_ch = t[index]
		bin_n, bin_edges = np.histogram(t_ch, bins=bins)
		bin_rate = bin_n / bg_dt
		bin_rate = np.concatenate((bin_rate[:1], bin_rate))
		bin_c, cs, bs = TD_baseline(bin_edges, bin_rate)
		if edges is not None:
			len_edges = edges[1:] - edges[:-1]
			bin_n1, bin_edges1 = np.histogram(t_ch, bins=edges)
			bin_c1 = (bin_edges1[1:] + bin_edges1[:-1]) * 0.5
			
		else:
			len_edges = bin_arr[:,-1]-bin_arr[:,0]
			bin_c1 = 0.5*(bin_arr[:,-1]+bin_arr[:,0])
			bin_n1 = overlap_curve(t_ch,bins_arr = bin_arr)
		bin_rate1 = bin_n1 / len_edges
		bs1 = np.interp(bin_c1, bin_c, bs)
		bin_err = np.sqrt(bin_n1) / len_edges
		bin_err[bin_err <= 0] = 1.0
		cs1 = bin_rate1 - bs1
		cs1[cs1 < 0] = 0
		data.append(list(cs1))
		data_err.append(list(bin_err))
	
	data = np.array(data).T
	data_err = np.array(data_err).T
	return data, data_err


