
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.patches import Polygon
from moviepy.editor import VideoClip
import numpy as np
from astropy.coordinates import SkyCoord,get_sun,get_body_barycentric,cartesian_to_spherical
from spherical_geometry.polygon import SphericalPolygon
import astropy.units as u


class Sky_map(object):

	def __init__(self,**kwargs):


		self.detector_color = '#74787c'
		self.detector_color_highlight = '#f26522'

		if 'figsize' in kwargs:
			self.size = kwargs['figsize'][0]
		else:
			kwargs['figsize'] = (5, 5)
			self.size = 5
		self.fig = plt.figure(**kwargs)
		self.ax = None



	def plot_earth(self, t, satellite, **kwargs):

		if 'facecolor' not in kwargs:
			kwargs['facecolor'] = '#90d7ec'
		if 'edgecolor' not in kwargs:
			kwargs['edgecolor'] = '#90d7ec'

		try:
			met, ra, dec, radius = satellite.get_earth_point(t)
			for i in range(len(met)):
				poly = SphericalPolygon.from_cone(ra[i],dec[i],radius[i],steps=180)
				x_,y_ = np.array(list(poly.to_radec())[0])
				earth = self.Polygon(list(zip(x_,y_))[::-1],**kwargs)
				self.ax.add_patch(earth)
		except (TypeError):

			met, ra, dec, radius = satellite.get_earth_point(t)
			poly = SphericalPolygon.from_cone(ra,dec,radius,steps=180)
			x_,y_ = np.array(list(poly.to_radec())[0])
			earth = self.Polygon(list(zip(x_,y_))[::-1],**kwargs)
			self.ax.add_patch(earth)


	def plot_detector(self, t,satellite,radius =10.0,good_detector_list = None,
			  detector_color = None,detector_color_highlight = None,
			  **kwargs):

		namelist = satellite.detectors.name
		try:
			len(t)
			print('only the first of t will be ploted')
			dete_point = satellite(t)
			for index,dete in enumerate(namelist):

				color = self.detector_color
				if detector_color is not None:
					try:
						color = detector_color[index]
					except (TypeError):
						color = detector_color

				if good_detector_list is not None:
					if dete in good_detector_list:
						color = self.detector_color_highlight
						if detector_color_highlight is not None:
							try:
								color = detector_color_highlight[index]
							except (TypeError):
								color = detector_color_highlight
				ra,dec = dete_point[dete][0]
				poly = SphericalPolygon.from_cone(ra,dec,radius,steps=180)
				x_,y_ = np.array(list(poly.to_radec())[0])
				dete_p = self.Polygon(list(zip(x_,y_))[::-1],facecolor = color,edgecolor=color,linewidth=0.2*self.size,alpha=0.5,**kwargs)
				self.ax.add_patch(dete_p)
				self.text(ra,dec,str(dete),color = 'k',va = 'center',ha='center',size=self.size,**kwargs)
		except (TypeError):
			dete_point = satellite(t)
			for index,dete in enumerate(namelist):
				color = self.detector_color
				if detector_color is not None:
					try:
						color = detector_color[index]
					except (TypeError):
						color = detector_color
				if good_detector_list is not None:
					if dete in good_detector_list:
						color = self.detector_color_highlight
						if detector_color_highlight is not None:
							try:
								color = detector_color_highlight[index]
							except (TypeError):
								color = detector_color_highlight
				ra,dec = dete_point[dete]
				poly = SphericalPolygon.from_cone(ra,dec,radius,steps=180)
				x_,y_ = np.array(list(poly.to_radec())[0])
				dete_p = self.Polygon(list(zip(x_,y_))[::-1],facecolor = color,edgecolor=color,linewidth=0.2*self.size,alpha=0.5,**kwargs)
				self.ax.add_patch(dete_p)
				self.text(ra,dec,str(dete),color = 'k',va = 'center',ha='center',size=self.size,**kwargs)


	def add_source(self, point, name=None):

		if isinstance(point,SkyCoord):
			ra = point.ra.deg
			dec = point.dec.deg
		else:
			ra = point[0]
			dec = point[1]
		self.plot(ra,dec,'*',color='r', markersize=self.size)
		if name is not None:
			try:
				for i in range(len(ra)):
					self.text(ra[i],dec[i],name[i],va = 'bottom',ha='right',size=self.size)
			except(TypeError):
				self.text(ra,dec,str(name),va = 'bottom',ha='right',size=self.size)

	def plot_continue_source(self):
		ra = [299.591,224.980,135.529,83.633]
		dec = [35.2020, -15.639, -40.5547, 22.0069]
		name = ['CygX-1','SocX-1','VeleX-1','Crab']
		self.plot(ra,dec,'o',color='#c7a252', markersize=self.size)
		for ind_,name_i in enumerate(name):
			self.text(ra[ind_],dec[ind_],name_i,va = 'bottom',ha='right',size=self.size)


	def plot_galactic_plane(self):

		l = np.linspace(0,360,720)
		b = np.zeros(l.size)
		sky = SkyCoord(l=l,b=b,frame = 'galactic',unit = 'deg')
		self.plot(sky.icrs.ra.deg,sky.icrs.dec.deg,color = '#281f1d',linewidth=0.4*self.size,alpha=0.3)
		self.plot(sky.icrs.ra.deg,sky.icrs.dec.deg,color = '#281f1d',linewidth=0.2*self.size,alpha=0.5)
		self.plot(sky.icrs.ra.deg[0],sky.icrs.dec.deg[0],'o',markersize = 0.4*self.size,color = '#281f1d',alpha=0.5)
		self.plot(sky.icrs.ra.deg[0],sky.icrs.dec.deg[0],'o',markersize = 0.6*self.size,color = '#281f1d',alpha=0.3)
		return sky.icrs.ra.deg,sky.icrs.dec.deg

	def plot_sum(self,t,satellite):

		clock = satellite.clock
		utc_t = clock.met_to_utc(t)
		tmp_sun = get_sun(utc_t)
		self.plot(tmp_sun.ra.deg,tmp_sun.dec.deg,'o',color='#ffd400',markersize=2*self.size)
		self.text(tmp_sun.ra.deg,tmp_sun.dec.deg,'sun',size=self.size,va = 'center',ha='center')
		return tmp_sun

	def plot_moon(self, t,satellite):
		clock = satellite.clock
		time_utc = clock.met_to_utc(t)
		earth_r = get_body_barycentric('earth', time_utc)
		moon_r = get_body_barycentric('moon',time_utc )
		r_e_m = moon_r - earth_r
		e = satellite.get_pos(t)
		r = -e*satellite.pos_unit - np.array([r_e_m.x.value, r_e_m.y.value, r_e_m.z.value]) * u.km
		moon_point_d = cartesian_to_spherical(-r[0], -r[1], -r[2])
		self.plot(moon_point_d[2].deg, moon_point_d[1].deg, 'o', color='#72777b', markersize=1.5*self.size)
		self.text(moon_point_d[2].deg, moon_point_d[1].deg, 'moon', size=self.size,va = 'center',ha='center')

	def title(self,*args,**kwargs):
		if self.ax is None:
			self.add_subplot(1,1,1)
		self.ax.set_title(*args,**kwargs)

	def add_subplot(self, *args, **kwargs):
		if 'lon_0' in kwargs:
			lon_0 = kwargs['lon_0']
		else:
			lon_0 = 180
		if 'facecolor' not in kwargs:
			kwargs['facecolor'] = '#f6f5ec'
		self.ax = self.fig.add_subplot(projection=ccrs.Mollweide(central_longitude=lon_0), *args, **kwargs)
		self.set_coordinates(lon_0, size=self.size)
		return self.ax

	def legend(self,*args, **kwargs):
		return self.ax.legend(*args, **kwargs)

	def plot(self,*args, **kwargs):
		if self.ax is None:
			self.add_subplot(1,1,1)
		return self.ax.plot(transform=ccrs.Geodetic(),*args, **kwargs)
	def scatter(self,x,y,*args, **kwargs):
		if self.ax is None:
			self.add_subplot(1,1,1)
		lats_y_ticke=self.ax.projection.transform_points(ccrs.Geodetic(),x, y)
		y_x = lats_y_ticke[:,0]
		y_y = lats_y_ticke[:,1]
		return self.ax.scatter(y_x,y_y,*args, **kwargs)
	def tricontour(self,*args, **kwargs):
		if self.ax is None:
			self.add_subplot(1,1,1)
		return self.ax.tricontour(transform=ccrs.Geodetic(),*args, **kwargs)
	def tricontourf(self,*args, **kwargs):
		if self.ax is None:
			self.add_subplot(1,1,1)
		return self.ax.tricontourf(transform=ccrs.Geodetic(),*args, **kwargs)
	def savefig(self,*args, **kwargs):
		return self.fig.savefig(*args, **kwargs)
	def close(self):
		plt.close(self.fig)
	def text(self,*args, **kwargs):
		if self.ax is None:
			self.add_subplot(1,1,1)
		return self.ax.text(transform=ccrs.Geodetic(),*args, **kwargs)
	def Polygon(self,*args, **kwargs):
		return Polygon(transform=ccrs.Geodetic(),*args, **kwargs)
	def add_patch(self,*args, **kwargs):
		return self.ax.add_patch(*args, **kwargs)

	def set_coordinates(self,lon_0,size = None):
		xticks = list(range(-180, 180, 30))
		yticks = list(range(-90, 90, 15))
		lons_x = np.arange(0,360,30)
		lons_y = np.zeros(lons_x.size)
		lats_y = np.arange(-75,76,15)
		lats_x = np.zeros(lats_y.size)
		self.ax.gridlines(xlocs=xticks, ylocs=yticks)
		lats_y_ticke = self.ax.projection.transform_points(ccrs.Geodetic(),lats_x+lon_0+179.99, lats_y*1.1)
		lats_y_x = lats_y_ticke[:,0]*0.86
		lats_y_y = lats_y_ticke[:,1]
		proj_xyz = self.ax.projection.transform_points(ccrs.Geodetic(),np.array([0,30]), np.array([0,0]))
		dx_ = np.abs(proj_xyz[0][0]-proj_xyz[1][0])
		for indexi,i in enumerate(lons_x):
			self.ax.text(i,lons_y[indexi],r'$%d^{\circ}$'%i,transform = ccrs.Geodetic(),size = size)
		for indexi,i in enumerate(lats_y):
			self.ax.text(lats_y_x[indexi]+dx_,lats_y_y[indexi],r'$%d^{\circ}$'%i,ha = 'right',va = 'center',size = size)
		self.ax.set_global()
		self.ax.invert_xaxis()


