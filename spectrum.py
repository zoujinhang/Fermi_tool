from astropy.io import fits
import numpy as np
from scipy import optimize
from Data_analysis.Baseline import WhittakerSmooth
#import json
import pymultinest
from .PlotMarginalModes import PlotMarginalModes
import matplotlib.pyplot as plt
import ctypes
import sys
paths = sys.path
c_lib_link = ''
for path in paths:
	fand = path+'/'

clib = ctypes.cdll.LoadLibrary('/home/laojin/my_lat/python_c/spectrum_tool.so')

class Fit(object):
	
	def __init__(self,spectrumlist,model,prior,parameters):
		'''
		
		:param spectrumlist:
		:param modellist:
		'''
		
		self.spec_num = len(spectrumlist)
		self.spectrumlist = spectrumlist
		self.model = model
		self.prior = prior
		self.parameters = parameters
		self.n_params = len(parameters)
		
		
		
	def log_like(self,cube,ndim,nparams):
		loglike = 0.
		#for spec_index in range(self.spec_num):
		for spec in self.spectrumlist:
			nn = len(spec.e_lo)
			spec1 = np.zeros(nn)
			#print('------------------------------------------------')
			#print('\r')
			for i in range(nn):
				e = np.linspace(spec.e_lo[i], spec.e_hi[i], 5)
				try:
					rate = self.model(e,cube)
					
					if True in np.isnan(rate):
						return -np.inf
					sp = rate.mean()
					spec1[i] = sp
				except:
					print('there are something wrong in your model!')
					print('the model`s return is ',cube)
					return -np.inf
			yu = spec.transform(spec1)
			if spec.effective_index is not None:
				loglike = loglike-0.5*(((yu[spec.effective_index[0]:spec.effective_index[-1]] -spec.spectrum[spec.effective_index[0]:spec.effective_index[-1]])/spec.spectrum_err[spec.effective_index[0]:spec.effective_index[-1]])**2).sum()
			else:
				loglike = loglike-0.5*(((yu - spec.spectrum)/spec.spectrum_err)**2).sum()
		return loglike
	def log_like_c(self,cube,ndim,nparams):
		loglike = 0.
		for spec in self.spectrumlist:
			e_add = spec.e_add
			try:
				rate = self.model(e_add,cube)
				if True in np.isnan(rate):
					return -np.inf
			except:
				print('there are something wrong in your model!')
				print('the model`s return is ',cube)
				return -np.inf
			spec1 = self.get_A(rate,spec.e_add_num)
			yu = spec.transform(spec1)
			if spec.effective_index is not None:
				loglike = loglike-0.5*(((yu[spec.effective_index[0]:spec.effective_index[-1]] -spec.spectrum[spec.effective_index[0]:spec.effective_index[-1]])/spec.spectrum_err[spec.effective_index[0]:spec.effective_index[-1]])**2).sum()
			else:
				loglike = loglike-0.5*(((yu - spec.spectrum)/spec.spectrum_err)**2).sum()
		return loglike

	def get_A(self,spe,add_n):
		n_spe = len(spe)
		spe = (ctypes.c_double * n_spe)(*list(spe))
		ret = (ctypes.c_double* int(n_spe/add_n))()
		clib.A_spec(spe,ret,n_spe,add_n,int(n_spe/add_n))
		return np.array(ret)
	
	def run(self,outputfiles_basename,resume = False, verbose = True):
		'''
		
		:param outputfiles_basename:
		:param resume:
		:param verbose:
		:return:
		'''
		pymultinest.run(self.log_like_c, self.prior, self.n_params, outputfiles_basename=outputfiles_basename,resume = resume, verbose = verbose)
		a1 = pymultinest.Analyzer(outputfiles_basename=outputfiles_basename, n_params = self.n_params)
		return a1
		
	def plot_model(self,a1,n = 0,ax = None,reference = True):
		'''
		
		:param a1:
		:param ax:
		:param reference:
		:return:
		'''
		e_c_list =[]
		for spec in self.spectrumlist:
			e_c = np.sqrt(spec.e_lo * spec.e_hi)
			e_c_list.append(e_c)
		for para in a1.get_equal_weighted_posterior()[::100,:-1]:
			for e_c in e_c_list:
				if ax is not None:
					ax.plot(e_c,e_c**n*self.model(e_c,para),'-', color='k', alpha=0.2)
				else:
					plt.plot(e_c,e_c**n*self.model(e_c,para),'-', color='k', alpha=0.2)
		best_value = a1.get_best_fit()['parameters']
		if reference:
			for index,spec in enumerate(self.spectrumlist):
				e_c,xxx = spec.get_reference()
				if ax is not None:
					ax.plot(e_c, e_c**n*xxx, '-.',label =spec.name +  ' reference')
				else:
					plt.plot(e_c,e_c**n*xxx , '-.',label =spec.name + ' reference')
				
		for e_c in e_c_list:
			if ax is not None:
				ax.plot(e_c,e_c**n*self.model(e_c,best_value) , '-',label = 'best model')
			else:
				plt.plot(e_c,e_c**n*self.model(e_c,best_value) , '-',label = 'best model')
				
		
	def plot_data(self,a1,ax = None,reference = True):
		'''
		
		:param a1:
		:param ax:
		:param reference:
		:return:
		'''
		best_value = a1.get_best_fit()['parameters']
		for spec in self.spectrumlist:
			e_c0 = spec.e_c
			e_c,xxx = spec.get_reference()
			sp = spec.spectrum
			sp_er = spec.spectrum_err
			effinde = spec.effective_index
			if effinde is not None:
				if ax is not None:
					ax.errorbar(e_c0[effinde[0]:effinde[-1]],sp[effinde[0]:effinde[-1]],yerr = sp_er[effinde[0]:effinde[-1]],elinewidth=2,capsize=2,label = spec.name + ' data',alpha=0.3)
					if reference:
						ax.plot(e_c0[effinde[0]:effinde[-1]],spec.transform(xxx)[effinde[0]:effinde[-1]], '-.',label = spec.name + ' reference')
					ax.plot(e_c0[effinde[0]:effinde[-1]],spec.transform(self.model(e_c,best_value))[effinde[0]:effinde[-1]],label = spec.name + ' model')
				else:
					plt.errorbar(e_c0[effinde[0]:effinde[-1]],sp[effinde[0]:effinde[-1]],yerr = sp_er[effinde[0]:effinde[-1]],elinewidth=2,capsize=2,label = spec.name + ' data',alpha=0.3)
					if reference:
						plt.plot(e_c0[effinde[0]:effinde[-1]],spec.transform(xxx)[effinde[0]:effinde[-1]], '-.',label = spec.name + ' reference')
					plt.plot(e_c0[effinde[0]:effinde[-1]],spec.transform(self.model(e_c,best_value))[effinde[0]:effinde[-1]],label = spec.name + ' model')
			else:
				if ax is not None:
					ax.errorbar(e_c0,sp,yerr = sp_er,elinewidth=2,capsize=2,label = spec.name + ' data',alpha=0.3)
					if reference:
						ax.plot(e_c0,spec.transform(xxx), '-.',label = spec.name + ' reference')
					ax.plot(e_c0,spec.transform(self.model(e_c,best_value)),label = spec.name + ' model')
				else:
					plt.errorbar(e_c0,sp,yerr = sp_er,elinewidth=2,capsize=2,label = spec.name + ' data',alpha=0.3)
					if reference:
						plt.plot(e_c0,spec.transform(xxx), '-.',label = spec.name + ' reference')
					plt.plot(e_c0,spec.transform(self.model(e_c,best_value)),label = spec.name + ' model')

			
	def plot_corner(self,a1):
		'''
		
		:param a1:
		:return:
		'''
		c = a1.get_stats()['modes'][0]
		mean = np.array(c['mean'])
		sigma = np.array(c['sigma'])
		maximum = np.array(c['maximum'])
		dx = mean-maximum
		sigmal = sigma-dx
		sigmah = sigma+dx
		p = PlotMarginalModes(a1)
		plt.figure(figsize=(5*self.n_params,5*self.n_params))
		for i in range(self.n_params):
			for j in range(i,self.n_params):
				if(i == 0):
					plt.subplot(self.n_params, self.n_params,i*self.n_params+j+1)
					plt.title(self.parameters[j],size = 30)
					plt.tick_params(labelsize = 15)
					p.plot_marginal(j, with_ellipses = True, with_points = False, grid_points=50)
					plt.errorbar(maximum[j],0,xerr=[[sigmal[j]],[sigmah[j]]],fmt = 'o',ecolor = 'r',capthick=2,color = 'r')
					if(j == 0):
						plt.ylabel("Probability",size = 30)
						plt.xlabel(self.parameters[j],size = 30)
					#plt.tight_layout()
				else:
					plt.subplot(self.n_params, self.n_params,i*self.n_params+j+1)
					plt.tick_params(labelsize = 15)
					p.plot_conditional(j, i-1, with_ellipses = False, with_points = False, grid_points=30)
					plt.errorbar(maximum[j], maximum[i-1], xerr=[[sigmal[j]],[sigmah[j]]], yerr=[[sigmal[i-1]],[sigmah[i-1]]], fmt='o', ecolor='r', capthick=2, color='r')
					if(j == i):
						plt.xlabel(self.parameters[j],size = 30)
						plt.ylabel(self.parameters[i-1],size = 30)
		
class Spectrum(object):
	'''
	
	'''
	def __init__(self,spectrum,spectrum_err,rsp_link,effective_band=None,spectrum_name = 'data',e_add_num = 10):
		hl = fits.open(rsp_link)
		matr = hl[2].data['MATRIX']
		matr[-1] = np.zeros(128)
		matr = np.vstack(matr)
		self.R = np.mat(matr)
		self.e_add_num = e_add_num
		self.name = spectrum_name
		self.e_lo = hl[2].data['ENERG_LO']#140
		self.e_hi = hl[2].data['ENERG_HI']
		self.e_add = self.get_e_add(self.e_lo,self.e_hi,self.e_add_num)
		self.e_min = hl[1].data['E_MIN']#128
		self.e_max = hl[1].data['E_MAX']
		self.e_c = np.sqrt(self.e_min*self.e_max)
		self.spectrum = np.array(spectrum)
		self.spectrum_err = np.array(spectrum_err)
		self.Y = np.mat(self.spectrum)
		self.Y_E = np.mat(self.spectrum_err)
		self.effective_band = effective_band
		if effective_band is not None:
			index_ = np.where((self.e_c>=effective_band[0])&(self.e_c<=effective_band[1]))[0]
			self.effective_index = [index_[0],index_[-1]+1]
		else:
			self.effective_index = None
	def get_e_add(self,e_lo,e_hi,nn = 5):
		e_add = np.array([])
		for i in range(len(e_lo)):
			e = np.linspace(e_lo[i], e_hi[i], nn)
			e_add = np.concatenate((e_add,e))
		return e_add
	def change_spectrum(self,spectrum,spectrum_err,effective_band=None):
		self.spectrum = np.array(spectrum)
		self.spectrum_err = np.array(spectrum_err)
		self.Y = np.mat(self.spectrum)
		self.Y_E = np.mat(self.spectrum_err)
		if effective_band is not None:
			index_ = np.where((self.e_c>=effective_band[0])&(self.e_c<=effective_band[1]))[0]
			self.effective_index = [index_[0],index_[-1]+1]
		
	def get_reference(self):
		x0 = np.zeros(len(self.e_lo))
		bounds = [[0,+np.inf]]*len(self.e_lo)
		xxx = optimize.fmin_l_bfgs_b(self.fmin,x0,fprime = self.fprime,bounds = bounds,epsilon=1.e-10)[0]
		xxx1 = WhittakerSmooth(xxx,x0+1,1)
		e_c = np.sqrt(self.e_hi*self.e_lo)
		return e_c,xxx1
		
	def transform(self,model):
		model_mat = np.mat(model)
		return np.array(model_mat*self.R)[0]
	
	
	def fprime(self,x):
		X = np.mat(x)
		G = 2*((X*self.R-self.Y)/self.Y_E/self.Y_E)*self.R.T
		return np.array(G)[0]
	
	
	def fmin(self,x):
		X = np.mat(x)
		M = (X*self.R-self.Y)/self.Y_E
		S = M*M.T
		return S[0,0]







