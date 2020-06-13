
import pandas as pd

class Save_track(object):
	
	def __init__(self,result,geometry):
		self.result = result
		self.geometry = geometry
		self.clock = geometry.Time_transition
		
	def save_bayesian_responses(self,sn,savename):
		sn_serch_result = self.result[sn]
		sn_trig_0 = sn_serch_result['trig_1']
		ls = sn_trig_0.shape[0]
		if ls > 0:
			sn_trig = sn_trig_0.drop_duplicates('start','first',inplace=False,ignore_index=True)
			start_met = sn_trig['start'].values
			stop_met = sn_trig['stop'].values
			start_utc = self.clock.batch_met_to_utc(start_met).fits
			stop_utc = self.clock.batch_met_to_utc(stop_met).fits
			overlap = sn_trig['overlap'].values
			c = {'start_utc':start_utc,
			     'stop_utc':stop_utc,
			     'start_met':start_met,
			     'stop_met':stop_met,
			     'overlap':overlap}
			hl = pd.DataFrame(c)
			hl.to_csv(savename, index=False)
		else:
			print(savename)
			print('the shape of bayesian_responses is 0.')
	
	def save_threshold_responses(self,sn,savename):
		sn_serch_result = self.result[sn]
		sn_trig_0 = sn_serch_result['trig_0']
		ls = sn_trig_0.shape[0]
		if ls > 0:
			sn_trig = sn_trig_0.drop_duplicates('start','first',inplace=False,ignore_index=True)
			start_met = sn_trig['start'].values
			stop_met = sn_trig['stop'].values
			start_utc = self.clock.batch_met_to_utc(start_met).fits
			stop_utc = self.clock.batch_met_to_utc(stop_met).fits
			overlap = sn_trig['overlap'].values
			c = {'start_utc':start_utc,
			     'stop_utc':stop_utc,
			     'start_met':start_met,
			     'stop_met':stop_met,
			     'overlap':overlap}
			hl = pd.DataFrame(c)
			hl.to_csv(savename, index=False)
		else:
			print(savename)
			print('The shape of threshold_responses is 0.')
	
	def save_all_respond(self,sn,savename):
		sn_serch_result = self.result[sn]
		sn_trig = sn_serch_result['trig_all']
		ls = sn_trig.shape[0]
		if ls > 0:
			start_met = sn_trig['start'].values
			stop_met = sn_trig['stop'].values
			start_utc = self.clock.batch_met_to_utc(start_met).fits
			stop_utc = self.clock.batch_met_to_utc(stop_met).fits
			detector = sn_trig['detector'].values
			bayes = sn_trig['bayes']
			c = {'start_utc':start_utc,
			     'stop_utc':stop_utc,
			     'start_met':start_met,
			     'stop_met':stop_met,
			     'detector':detector,
			     'bayes':bayes}
			hl = pd.DataFrame(c)
			hl.to_csv(savename, index=False)
		else:
			print(savename)
			print('The shape of all_respond is 0.')
			
			
class Save_search(object):
	
	def __init__(self,result,geometry):
		
		self.result = result
		self.geometry = geometry
		self.clock = geometry.Time_transition
		
	def save_bayesian_responses(self,savename):
		sn_trig_0 = self.result['trig_1']
		ls = sn_trig_0.shape[0]
		if ls > 0:
			sn_trig = sn_trig_0.drop_duplicates('start','first',inplace=False,ignore_index=True)
			start_met = sn_trig['start'].values
			stop_met = sn_trig['stop'].values
			start_utc = self.clock.batch_met_to_utc(start_met).fits
			stop_utc = self.clock.batch_met_to_utc(stop_met).fits
			overlap = sn_trig['overlap'].values
			c = {'start_utc':start_utc,
			     'stop_utc':stop_utc,
			     'start_met':start_met,
			     'stop_met':stop_met,
			     'overlap':overlap}
			hl = pd.DataFrame(c)
			hl.to_csv(savename, index=False)
		else:
			print(savename)
			print('the shape of bayesian_responses is 0.')
	
	def save_threshold_responses(self,savename):
		
		sn_trig_0 = self.result['trig_0']
		ls = sn_trig_0.shape[0]
		if ls > 0:
			sn_trig = sn_trig_0.drop_duplicates('start','first',inplace=False,ignore_index=True)
			start_met = sn_trig['start'].values
			stop_met = sn_trig['stop'].values
			start_utc = self.clock.batch_met_to_utc(start_met).fits
			stop_utc = self.clock.batch_met_to_utc(stop_met).fits
			overlap = sn_trig['overlap'].values
			c = {'start_utc':start_utc,
			     'stop_utc':stop_utc,
			     'start_met':start_met,
			     'stop_met':stop_met,
			     'overlap':overlap}
			hl = pd.DataFrame(c)
			hl.to_csv(savename, index=False)
		else:
			print(savename)
			print('The shape of threshold_responses is 0.')
	
	def save_all_respond(self,savename):
		
		sn_trig =  self.result['trig_all']
		ls = sn_trig.shape[0]
		if ls > 0:
			start_met = sn_trig['start'].values
			stop_met = sn_trig['stop'].values
			start_utc = self.clock.batch_met_to_utc(start_met).fits
			stop_utc = self.clock.batch_met_to_utc(stop_met).fits
			detector = sn_trig['detector'].values
			bayes = sn_trig['bayes']
			c = {'start_utc':start_utc,
			     'stop_utc':stop_utc,
			     'start_met':start_met,
			     'stop_met':stop_met,
			     'detector':detector,
			     'bayes':bayes}
			hl = pd.DataFrame(c)
			hl.to_csv(savename, index=False)
		else:
			print(savename)
			print('The shape of all_respond is 0.')





