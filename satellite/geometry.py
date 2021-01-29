
from .detectors import Detectors
from .clock import Clock
import numpy as np
import astropy.units as u
from scipy.interpolate import interp1d
from astropy.coordinates import cartesian_to_spherical,SkyCoord
import pandas as pd

class Geometry(object):

	def __init__(self,pd_position_data,pos_unit = u.m,detectors = None, clock = None):

		if detectors is not None:
			self.detectors = detectors
		else:
			self.detectors = Detectors()

		if clock is not None:

			self.clock = clock
		else:
			self.clock = Clock()

		#self.pd_position_data = pd_position_data
		met_time = pd_position_data['SCLK_UTC'].values
		self.met_time_band = [met_time.min(), met_time.max()]
		met_time[0] = met_time[0] - 0.5
		met_time[-1] = met_time[-1] + 0.5


		qsj = pd_position_data[['QSJ_1', 'QSJ_2', 'QSJ_3', 'QSJ_4']].values
		self.deter_f = {}
		mat_list = self.get_transform_mat_more(qsj)
		for detec_i in self.detectors.get_names():
			vector = self.detectors(detec_i)
			position_list = []
			for mat in mat_list:
				X = np.mat(vector).T
				X1 = np.array(mat * X).T[0]
				position_list.append([X1[0],X1[1],X1[2]])
			position_list = np.array(position_list).T
			x_f = interp1d(met_time,position_list[0],kind = 'quadratic')
			y_f = interp1d(met_time, position_list[1], kind='quadratic')
			z_f = interp1d(met_time, position_list[2], kind='quadratic')
			self.deter_f[detec_i] = [x_f,y_f,z_f]
		self.pos_unit = pos_unit
		self.qsj0_f = interp1d(met_time,pd_position_data['QSJ_1'].values,kind='quadratic')
		self.qsj1_f = interp1d(met_time,pd_position_data['QSJ_2'].values,kind='quadratic')
		self.qsj2_f = interp1d(met_time,pd_position_data['QSJ_3'].values,kind='quadratic')
		self.qsj3_f = interp1d(met_time,pd_position_data['QSJ_4'].values,kind='quadratic')
		self.pos_x_f = interp1d(met_time,-pd_position_data['POS_X'].values,kind='quadratic')
		self.pos_y_f = interp1d(met_time,-pd_position_data['POS_Y'].values,kind='quadratic')
		self.pos_z_f = interp1d(met_time,-pd_position_data['POS_Z'].values,kind='quadratic')

	def __call__(self, met_time):

		try:
			len(met_time)
			met_time = np.array(met_time)
			index_in = np.where((met_time >= self.met_time_band[0]) & (met_time <= self.met_time_band[-1]))[0]
			met_time_in = met_time[index_in]
			c = {'time': met_time_in}
			for detec in self.detectors.get_names():
				xf, yf, zf = self.deter_f[detec]
				x = xf(met_time_in)
				y = yf(met_time_in)
				z = zf(met_time_in)
				position = cartesian_to_spherical(x, y, z)
				c[detec] = np.vstack([position[2].deg, position[1].deg]).T
			return c

		except (TypeError):

			c = {'time': met_time}
			for detec in self.detectors.get_names():
				xf, yf, zf = self.deter_f[detec]
				x = xf(met_time)
				y = yf(met_time)
				z = zf(met_time)
				position = cartesian_to_spherical(x, y, z)
				c[detec] = [position[2].deg, position[1].deg]
			return c

	def get_separation_with_time(self, time, source):

		name = self.detectors.get_names()
		try:
			len(time)
			time = np.array(time)
			index_in = np.where((time >= self.met_time_band[0]) & (time <= self.met_time_band[-1]))[0]
			time_in = time[index_in]
			dete_centor = self(time_in)
			data = [time_in]
			col_name = ['time']
			for ni in name:
				center = (dete_centor[ni]).T
				center = SkyCoord(ra=center[0], dec=center[1], frame='icrs', unit='deg')
				data.append(center.separation(source).value)
				col_name.append(ni)
			data = np.vstack(data).T
		except(TypeError):
			data = [time]
			col_name = ['time']
			dete_centor = self(time)
			for ni in name:
				center = (dete_centor[ni])
				center = SkyCoord(ra=center[0], dec=center[1], frame='icrs', unit='deg')
				data.append(center.separation(source).value)
				col_name.append(ni)
			data = np.array([data])
		return pd.DataFrame(data=data, columns=col_name, dtype=np.float128)


	def get_earth_point(self,met_time):

		earth_radius = 6371. * u.km
		try:
			len(met_time)
			met_time = np.array(met_time)
			index_in = np.where((met_time >= self.met_time_band[0]) & (met_time <= self.met_time_band[-1]))[0]
			met_time_in = met_time[index_in]
			pos_x = self.pos_x_f(met_time_in)* self.pos_unit
			pos_y = self.pos_y_f(met_time_in)* self.pos_unit
			pos_z = self.pos_z_f(met_time_in)* self.pos_unit
			position = cartesian_to_spherical(pos_x,pos_y,pos_z)
			fermi_radius = np.sqrt(pos_x**2 + pos_y**2 + pos_z**2)
			radius_deg = np.rad2deg(np.arcsin((earth_radius / fermi_radius).to(u.dimensionless_unscaled)).value)
			return met_time_in,position[2].deg,position[1].deg,radius_deg

		except (TypeError):

			pos_x = self.pos_x_f(met_time)* self.pos_unit
			pos_y = self.pos_y_f(met_time)* self.pos_unit
			pos_z = self.pos_z_f(met_time)* self.pos_unit
			position = cartesian_to_spherical(pos_x,pos_y,pos_z)
			fermi_radius = np.sqrt(pos_x**2 + pos_y**2 + pos_z**2)
			radius_deg = np.rad2deg(np.arcsin((earth_radius / fermi_radius).to(u.dimensionless_unscaled)).value)
			return met_time,position[2].deg,position[1].deg,radius_deg


	def get_transform_mat_more(self,qsj):

		mat_list = []
		for qsj_i in qsj:
			mat_list.append(self.get_transform_mat_one(qsj_i))
		return mat_list

	def get_pos(self,time):

		try:
			len(time)
			pos_x = self.pos_x_f(time)
			pos_y = self.pos_y_f(time)
			pos_z = self.pos_z_f(time)
			pos = np.vstack([pos_x,pos_y,pos_z])
			return pos.T

		except (TypeError):

			return np.array([self.pos_x_f(time),self.pos_y_f(time),self.pos_z_f(time)])

	def get_qsj(self,time):

		try:
			len(time)
			qsj0 = self.qsj0_f(time)
			qsj1 = self.qsj1_f(time)
			qsj2 = self.qsj2_f(time)
			qsj3 = self.qsj3_f(time)
			qsj = np.vstack([qsj0,qsj1,qsj2,qsj3])
			return qsj.T

		except (TypeError):

			return np.array([self.qsj0_f(time),self.qsj1_f(time),self.qsj2_f(time),self.qsj3_f(time)])


	def get_transform_mat_one(self,qsj):

		p1,p2,p3,p0 = qsj
		mat = np.mat(np.zeros((3,3)))
		mat[0, 0] = p0 ** 2 + p1 ** 2 - p2 ** 2 - p3 ** 2
		mat[0, 1] = 2 * (p1 * p2 - p0 * p3)
		mat[0, 2] = 2 * (p0 * p2 + p1 * p3)
		mat[1, 0] = 2 * (p3 * p0 + p2 * p1)
		mat[1, 1] = p0 ** 2 + p2 ** 2 - p3 ** 2 - p1 ** 2
		mat[1, 2] = 2 * (p2 * p3 - p1 * p0)
		mat[2, 0] = 2 * (p1 * p3 - p0 * p2)
		mat[2, 1] = 2 * (p0 * p1 + p3 * p2)
		mat[2, 2] = p0 ** 2 + p3 ** 2 - p1 ** 2 - p2 ** 2
		return mat
	def transform_frame_more(self,vector3d,qsj):
		matr = self.get_transform_mat_one(qsj)
		X = np.mat(vector3d).T
		return np.array(matr * X).T

	def transform_frame_one(self,vector3d,qsj):
		matr = self.get_transform_mat_one(qsj)
		X = np.mat(vector3d).T
		return np.array(matr * X).T[0]

	def inv_transform_frame_one(self,vector3d,qsj):
		matr = self.get_transform_mat_one(qsj)
		X = np.mat(vector3d).T
		return np.array(matr.I * X).T[0]


	def get_detector_centers(self,met_time):

		try:
			len(met_time)
			met_time = np.array(met_time)
			index_in = np.where((met_time>=self.met_time_band[0])&(met_time<=self.met_time_band[-1]))[0]
			met_time_in = met_time[index_in]
			qsj = self.get_qsj(met_time_in)

			c = {'time':met_time_in}

			for detec in self.detectors.get_names():
				c[detec] = []
				local_dete_ver = self.detectors(detec)
				for qsj_i in qsj:
					dete_ver = self.transform_frame_one(local_dete_ver,qsj_i)
					position = cartesian_to_spherical(dete_ver[0],dete_ver[1],dete_ver[2])
					c[detec].append([position[2].deg,position[1].deg])
				c[detec] = np.array(c[detec])

			return c

		except (TypeError):

			qsj = self.get_qsj(met_time)
			c = {'time':met_time}
			for detec in self.detectors.get_names():
				local_dete_ver = self.detectors(detec)
				dete_ver = self.transform_frame_one(local_dete_ver,qsj)
				position = cartesian_to_spherical(dete_ver[0],dete_ver[1],dete_ver[2])
				c[detec] = [position[2].deg,position[1].deg]
			return c









