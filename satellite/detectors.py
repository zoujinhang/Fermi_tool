

import numpy as np
from astropy.coordinates import spherical_to_cartesian,SkyCoord


class Detectors(object):

	def __init__(self,local_az=None,local_zen=None,namelist = None,eff_angle=None):

		if local_az is None or local_zen is None or namelist is None or eff_angle is None:

			self.name = np.array(['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb','b0','b1'])
			self.detetor_vector = np.array(
				[[0.2446677589, 0.2523893824, 0.9361823057],
				 [0.5017318971, 0.5036621127, 0.7032706462],
				 [0.5233876659, 0.8520868147, -0.0036651682],
				 [0.5009495177, -0.5032279093, 0.7041386753],
				 [0.5468267487, -0.8372325378, -0.0047123847],
				 [0.9982910766, 0.0584352143, 0.0005236008],
				 [-0.2471260191, -0.2465229020, 0.9370993606],
				 [-0.5135631636, -0.5067957667, 0.6923950822],
				 [-0.5503349679, -0.8349438131, 0.0005235846],
				 [-0.5064476761, 0.5030998708, 0.7002865795],
				 [-0.5552650628, 0.8316411478, -0.0073303046],
				 [-0.9978547710, -0.0652279514, -0.0055850266],
				 [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]
			)
			self.eff_angle = np.array([60,60,60,60,60,60,60,60,60,60,60,60,90,90])

		else:
			self.name = np.array(namelist)
			self.eff_angle = np.array(eff_angle)
			self.detetor_vector = np.array(spherical_to_cartesian(1,local_az,local_zen)).T
		self.overlap = self.get_angle_overlap()


	def __call__(self,name):

		n_index = np.where(self.name == name)[0]
		detector_vector = self.detetor_vector[n_index]

		return detector_vector[0]

	def get_eff_angle(self,name):

		n_index = np.where(self.name == name)[0]
		eff_angle = self.eff_angle[n_index]
		return eff_angle[0]

	def get_names(self):
		return self.name

	def get_angle_overlap(self):
		position_array = self.detetor_vector.T
		center_all = SkyCoord(x=position_array[0],y=position_array[1],z=position_array[2],frame = 'icrs',representation='cartesian')
		c = {}
		for index,di in enumerate(self.name):
			centeri = center_all[index]
			seq = center_all.separation(centeri).deg
			dindx = (seq>0)&(seq- 1.*self.eff_angle<=10)
			c[di] = self.name[dindx]
		return c


	def detector_association(self,namelist,n = 3,m = 3):

		v_ni_list = []
		for ni in namelist:
			over_list = list(self.overlap[ni])
			over_list.append(ni)
			v_ni_list = v_ni_list+over_list
		ni_n = list(set(v_ni_list))
		ni_n = np.array(ni_n)
		ni_nm =[]
		v_ni_list = np.array(v_ni_list)
		for ni in ni_n:
			index_n = np.where(v_ni_list == ni)[0]
			nn = v_ni_list[index_n].size
			ni_nm.append(nn)
		ni_nm = np.array(ni_nm)
		index_ = np.where(ni_nm>=n)[0]
		ni_n = ni_n[index_]
		ni_over_set = set(ni_n)
		ni_set = set(namelist)
		ni_union_set = list(ni_set & ni_over_set)
		num = len(ni_union_set)
		return num >= m

	def select_points(self,point,detector_list,n = 1):

		num = np.zeros(len(point),dtype = int)
		if len(detector_list)<n:
			n = len(detector_list)
		for dete in detector_list:
			v = self(dete)
			#print('v',dete,v)
			v_xyz = SkyCoord(x=v[0],y=v[1],z=v[2],frame='icrs',representation='cartesian')
			#print('v_xyz',dete,v_xyz)
			#print('v_xyz eff',dete,self.get_eff_angle(dete))
			seq = point.separation(v_xyz).deg <=self.get_eff_angle(dete)+10
			num[seq] = num[seq]+1
		return np.where(num>=n)[0]