from astropy.time import Time

class Clock(object):

	def __init__(self,time_origin = None):

		if time_origin is None:
			self.time_origin_gps = 662342413.0
		elif isinstance(time_origin, str):
			self.time_origin_gps = Time(time_origin).gps
		else:
			self.time_origin_gps = Time(time_origin,format = 'mjd',scale = 'utc').gps

	def utc_to_met(self, utc, format=None, scale=None):

		if isinstance(utc, Time):
			gps = utc.gps
		else:
			try:
				gps = Time(utc, format=format, scale=scale).gps
			except (ValueError):
				gps = Time(utc, format='mjd', scale=scale).gps
		return gps - self.time_origin_gps

	def met_to_utc(self, met):

		gps = self.time_origin_gps + met

		return Time(gps, format='gps', scale='utc')





