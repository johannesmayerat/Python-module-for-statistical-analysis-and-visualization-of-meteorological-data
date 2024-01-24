# Shared python functions.
# 
# AUTHOR: Johannes Mayer (johannes.mayer@univie.ac.at, johannes.mayer.private@gmail.com)
# 
# VERSION: 2023-11-07
# 
# FUNCTIONS:
#	_rcparams
#	EffectiveSampleSize
#	Compute_Trend
#	Compute_SeasonalMean
#

_currentroutine = 'shared.py'

def _rcparams():
	import matplotlib
	
	fs_all = 27
	fw_all = 'normal' # light, normal
	
	matplotlib.rcParams['font.size'] = fs_all
	matplotlib.rcParams['font.weight'] = fw_all
	matplotlib.rcParams['legend.fontsize'] = fs_all 
	matplotlib.rcParams['legend.fancybox'] = False
	matplotlib.rcParams['legend.frameon'] = False
	matplotlib.rcParams['axes.labelweight'] = fw_all
	matplotlib.rcParams['axes.titleweight'] = fw_all
	matplotlib.rc('font', family='sans-serif')#, serif='cm10')
	matplotlib.rc('text', usetex=False)
	#matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
	#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
	
	return fs_all, fw_all
	
def EffectiveSampleSize(sample=[],ac='ar1',CL=0.95):
	"""
	This function computes the effective sample size based on the autocorrelation function of a timeseries.
	
	Input:
		sample (list) input sample for which the effective sample size is computed
		ac (str) type of computation, possible arguments are 'ar1' and 'full'
		CL (float) Confidence interval, default is 0.95
	
	Notes:
		* The computation with ac='full' is based on Loeb et al. (2022), DOI:10.1029/2022JD036686
		* The computation with ac='ar1' considers only the autocorrelation function at lag=1
	"""

	import numpy as np
	import scipy.stats as stats

	if len(sample) == 0:
		print(EffectiveSampleSize.__doc__)
		return

	if ac not in ['ar1','full']:
		print("Wrong argument for ac. Possible arguments are 'ar1' and 'full'")
		return

	Ntrue = len(sample)
	if ac == 'ar1':
		alpha = np.corrcoef(sample[1:],sample[:-1])[0,1]
		Neff = Ntrue*(1-alpha)/(1+alpha)
		return Neff
		
	elif ac == 'full':
		alpha = []
		for lag in range(1,Ntrue-1):
			alpha += [np.corrcoef(sample[lag:],sample[:-lag])[0,1]]
		
		alpha = np.array(alpha)
		
		sigma_k_squared = [1/Ntrue]
		for i in range(1,len(alpha)):
			sigma_k_squared += [1/Ntrue*(1+2*np.sum(alpha[:i]**2.))]
			
		sigma_k_squared = np.array(sigma_k_squared)
		
		
		alpha_cl = 1 - CL
		tval = stats.t.ppf(1. - alpha_cl/2.,Ntrue-1)
		
		alpha_signif = np.copy(alpha)
		sel_m = 0
		for m in range(len(alpha)-2):
			if alpha[m+1] < 0 and alpha[m+1] + alpha[m+2] < 0:
				sel_m = m + 1
				break
		
		for m in range(sel_m):
			if tval*np.sqrt(sigma_k_squared[m]) > alpha[m] > -tval*np.sqrt(sigma_k_squared[m]):
				alpha[m] = 0.
		
		Neff = Ntrue/(1+2*np.sum(alpha[:sel_m]))

		return Neff

		
def Compute_Trend(DateArr=[],idata=[],factor=120):
	"""
	This routine computes the linear trend of an 1D, 2D, or, 3D input field.
	
	Input: 
		DateArr (list of datetime objects)
		idata (array of floats) input data, can have all dimenions <= 3
		factor (int) trend factor which is multiplied with the trend coefficient, default is 120 months (1 decade)
		
	Output: 
		trend*factor, trend_intercept
		
	Notes:
		* This function can be used for gridded data as well as timeseries.
	"""
	
	from sklearn import linear_model
	import numpy as np
	
	if len(DateArr) == 0 or len(idata) == 0:
		print(Compute_Trend.__doc__)
		return
	
	ndims = idata.ndim
	shape = idata.shape
	model_ols =  linear_model.LinearRegression() 
	time = np.arange(len(DateArr)).reshape(-1,1)

	if ndims == 3 and shape[2] == len(DateArr):
		mask = idata.mask[:,:,0]		# wähle maske basierend auf 1. zeitschritt aus
		index = np.where(mask == False)	# extrahiere alle indizes, welche nicht maskiert sind,
		idata = idata[index].T			# reduziere daten auf ausgewählte indizes. letzte dimension (zeit) bleibt erhalten
		
	elif ndims == 2 and shape[1] == len(DateArr):
		idata = idata.T
		
	# ndims == 1 -> not necessary, see below

	model_ols.fit(time, idata) 

	if ndims == 3:
		trend_out = np.ma.zeros([shape[0],shape[1]])	# erstelle globales gitter, aber mit maske
		trend_out.mask = [True]							# notwendig weil np.ma.zeros kein 'mask=' argument hat
		trend_out.mask[index] = False					# entferne maske von allen gitterpunkten, die verwendet werden
		trend_out[index] = model_ols.coef_[:,0]*factor	# füge daten in output array ein
    	
		interc_out = np.ma.zeros([shape[0],shape[1]])
		interc_out.mask = [True]
		interc_out.mask[index] = False
		interc_out[index] = model_ols.intercept_*factor    	
		
		return trend_out, interc_out
		
	elif ndims == 2: 
		return model_ols.coef_.reshape([shape[0]])*factor, model_ols.intercept_.reshape([shape[0]])
		
	elif ndims == 1: 
		return model_ols.coef_[0]*factor, model_ols.intercept_		
	
	
	
def Compute_SeasonalMean(DateArr=[],idata=[],season='DJF'):
	"""
	This function computes the seasonal mean of a timeseries (1D) or gridded data (3D).
	
	Input:
		DateArr (list of datetime objects) time
		idata (list or array) input data, can be a timeseries or gridded data
		season (str) default is 'DJF'
		
	Output:
		DateArr_season, idata_season
	
	Notes: 
		* DateArr must start with central_month!
	"""

	import pandas as pd
	import numpy as np
	import matplotlib.pylab as plt

	if len(DateArr) == 0 or len(idata) == 0:
		print(Compute_SeasonalMean.__doc__)
		return

	Months = 'NDJFMAMJJASOND'
	DaysPerMonth = np.array([30,31,31,28,31,30,31,30,31,31,30,31,30,31])
	IndexMonth = np.array([11,12,1,2,3,4,5,6,7,8,9,10,11,12])
		
	nmonth = len(season)
	nh = int((nmonth-1)//2)
		
	idata = np.ma.array(idata)
	
	if idata.ndim == 1:
		idata = idata[np.newaxis,np.newaxis,:]
		
	time_dim = idata.shape[2]		
		
	if len(DateArr) != time_dim:
		print('ERROR! Inconsistent length.')
		return False	
		
	if nmonth % 2 == 0:
		print("ERROR! Wrong length of season. The argument 'season' must contain an odd number of months.")
		return False	

	try:
		DPM_season = DaysPerMonth[Months.index(season):Months.index(season)+nmonth] 
		central_month = IndexMonth[Months.index(season)+nh]
	except:
		print('ERROR! Undefined season.')
		return False


	DateArr_season = []
	for i in range(time_dim):
		if DateArr[i].month == central_month:
			DateArr_season += [DateArr[i]]		

	odata = np.ma.zeros([idata.shape[0],idata.shape[1],len(DateArr_season)])
	cnto = 0
	
	for i in range(time_dim):
		if 'F' in season:
			if DateArr[i].year%4 == 0:
				DPM_season[season.index('F')] = 29
			else:
				DPM_season[season.index('F')] = 28
		    
		if DateArr[i].month == central_month:	
			if i == 0:
				odata[:,:,cnto] = np.sum(idata[:,:,i:i+nh+1]*DPM_season[nh:],axis=2)/np.sum(DPM_season[nh:])
			elif i == time_dim - 1:
				odata[:,:,cnto] = np.sum(idata[:,:,i-nh:i+1]*DPM_season[:nh+1],axis=2)/np.sum(DPM_season[:nh+1])
			else:
				odata[:,:,cnto] = np.sum(idata[:,:,i-nh:i+nh+1]*DPM_season[:],axis=2)/np.sum(DPM_season)
				
			cnto += 1 
			
	odata = np.squeeze(odata)

	return DateArr_season, odata	
	
	
	
	
	
		
	
	
	

