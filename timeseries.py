# Python functions for timeseries.
# 
# AUTHOR: Johannes Mayer (johannes.mayer@univie.ac.at, johannes.mayer.private@gmail.com)
# 
# VERSION: 2023-11-07
# 
# FUNCTIONS:
#	ReadData
#	SaveData
#	PlotData
#	Apply_SNHT
#	Apply_MovingAverage
#	Compute_Tendency
#	Compute_Correlation
#	Compute_AnomalyClimatology
#	Compute_ConfidenceIntervals
#

_currentroutine = 'timeseries.py'
_default_opath_plot = './'
_default_opath_data = './'

from shared import *	
	
def ReadData(ifile=None,fmt='%Y-%m-%d'):
	"""
	Input: 
		ifile (str) file path, 
		fmt   (str) datetime format
	
	Output:
		datetime object, data, header 
		
	Notes:
		* File must contain 2 columns separated by one space
		* First column contains the date in the format fmt (e.g.,YYYY-MM-DD)
		* Second column contains the data
	"""
	
	import numpy as np
	import datetime

	if ifile == None:
		print(f'In routine: {_currentroutine}')
		print(ReadData.__doc__)
		return

	iDateArr = []
	DateArr = []
	header = []
	f = open(ifile,"r") 
	lines = f.readlines()

	hlen = 0
	while '#' in lines[hlen]:
		header += [float(lines[hlen].replace('#',''))]
		hlen += 1

	if len(header) == 0: header = [None]

	idata = np.ma.zeros(len(lines)-hlen)
	idata.mask = [False] # creates mask with False for all values -> maske same dimension as data
	for i in range(len(idata)):
		linestrings = lines[i+hlen].split(' ')
		iDateArr += [linestrings[0]]
		if linestrings[1].strip() == '--':
			idata.mask[i] = True
		else:
			idata[i] = float(linestrings[1])
	f.close()

	#for i in range(len(iDateArr)):
	#	DateArr += [datetime.datetime.strptime(iDateArr[i],fmt)]
		
	DateArr = [datetime.datetime.strptime(iDateArr[i],fmt) for i in range(len(iDateArr))]

	return DateArr, idata, np.array(header)
	

def SaveData(DateArr=[],idata=[],ofile=[],opath=_default_opath_data, dfmt='%Y-%m-%d',header=None):
	"""
	This function saves a timeseries to a .dat-file.
	
	Input:
		idata (list/array of floats) input data
		DateArr (list of datetime objects) time
		ofile (str) file name
		opath (str) defines output folder, default given by _default_opath_data
	"""
	
	if len(DateArr) == 0 or len(idata) == 0 or len(ofile) == 0:
		print(f'In routine: {_currentroutine}')
		print(SaveData.__doc__)
		return

	DateArr_fct = DateArr.copy()

	### FOR MONTHLY DATA, CHECK IF DAY = 15 
	if len(dfmt) <= 8 and DateArr_fct[0].day != 15 and DateArr_fct[0].hour == 0 and DateArr_fct[0].minute == 0:	
		for i in range(len(DateArr_fct)):
			DateArr_fct[i] = DateArr_fct[i].replace(day=15)
		
	### WRITE OUTPUT
	ofile = opath + ofile
	print(f' -- Save timeseries: {ofile}')
	file1 = open(ofile,"w") 
	if header != None:
		for i in range(len(header)):
			file1.write(f'# {header[i]}\n')
	for i in range(len(idata)):
		file1.write(f'{DateArr_fct[i].strftime(dfmt)} {idata[i]}\n')
	file1.close() 	
	
	

def PlotData(DateArr=[],idata=[], ilabel=None,\
		pfmt='%Y', title=None, figsize=(18,4), units=None, ncol=1, \
		icolor=['b','r','k','g','y'], ils=['-']*5, ilw=[2]*5, xlim=[], ylim=None,\
		trend=False, trend_multiplier=1,\
		savefig=False, opath=_default_opath_plot):
		
	"""
	Input:
		DateArr (list of datetime objects lists)
		idata (list of lists) input data
		
	Optional input: 
		ilabel (list of str) labels of each input
		pfmt (str) plot format of the date
		savefig (str) output filename
		...
	"""
	
	import numpy as np
	import scipy.stats as stats
	import matplotlib.pylab as plt
	import matplotlib.dates as mdates
	import datetime
	from sklearn import linear_model
	import warnings
	warnings.filterwarnings("ignore",category=UserWarning)
	
	from shared import _rcparams
	fs_all, fw_all = _rcparams()
	
	myFmt = mdates.DateFormatter(pfmt)

	if trend: model_ols =  linear_model.LinearRegression() 

	if len(DateArr) == 0 or len(idata) == 0:
		print(f'In routine: {_currentroutine}')
		print(PlotData.__doc__)
		return

	if ilabel == None:
		dlabel = ['-']*len(idata)
	else:
		dlabel = ilabel

	if xlim == []:
		xlim = [DateArr[0][0],DateArr[0][-1]]
	elif isinstance(xlim[0], str) == True:
		xlim[0] = datetime.datetime.strptime(xlim[0],'%Y-%m-%d')
		xlim[1] = datetime.datetime.strptime(xlim[1],'%Y-%m-%d')
		
	fig = plt.figure(figsize=figsize)
	ax = plt.subplot()
	
	ax.xaxis.set_major_formatter(myFmt)
	if title != None: plt.title(title)
	
	if units != None: 
		PlotUnits = ''
		for i in units.split(' '): # 
			estr = ' '
			if '**' in i: estr = '}$ '		
			PlotUnits += i.replace('**',r'$^{')+estr

		plt.ylabel(PlotUnits)
	
	if trend:
		SlopeArr = []
		
	for i in range(len(idata)):
		plt.plot(DateArr[i],idata[i],label=dlabel[i],color=icolor[i],ls=ils[i], lw=ilw[i], alpha=1.0)

		if trend:
			time_numeric = np.arange(len(DateArr[i])).reshape(-1,1)

			model_ols.fit(time_numeric, idata[i]) 
			slope = model_ols.coef_[0]
        		
			SlopeArr += [slope]
			y_trend = model_ols.intercept_ + slope*time_numeric
			plt.plot(DateArr[i],y_trend,color=icolor[i],ls='--', lw=ilw[i],alpha=0.5)
			
	if ilabel != None:		
		leg = fig.legend(loc='lower left',ncol=ncol)#,bbox_to_anchor=(1.0,1.25))
		plt.draw() # Draw the figure so you can find the positon of the legend. 
		bb = leg.get_bbox_to_anchor().inverse_transformed(ax.transAxes) # Get the bounding box of the original legend

		# Change to location of the legend. 
		yoffset = 1.10
		xoffset = 0.15
		bb.y0 += yoffset
		bb.y1 += yoffset
		bb.x0 += xoffset
		bb.x1 += xoffset
		leg.set_bbox_to_anchor(bb, transform=ax.transAxes)
	
	plt.xlim(xlim[0],xlim[1])
	if ylim != None: plt.ylim(ylim)
	plt.grid(ls=':')
	if savefig != False:
		ofile = opath+savefig
		plt.savefig(ofile,dpi=300,bbox_inches='tight',transparent=False,facecolor='white')
		print(f' --savefig: {ofile}')
	plt.show()
	
	if trend:
		for i in range(len(idata)):
			print(f'Trend {dlabel[i]}: {SlopeArr[i]*trend_multiplier}')
		return np.array(np.array(SlopeArr)*trend_multiplier)
	else: 
		return


def Apply_SNHT(idata=[], window=48):
	"""
	The function applies a standard normal homogeneity test on the input data. 
	
	Input: 
		idata (list of floats)
		window (int) SNHT window on one side, default is 48
		
	Output: 
		Tk (list of floats)
	"""
	
	import numpy as np

	if len(idata) == 0:
		print(f'In routine: {_currentroutine}')
		print(Apply_SNHT.__doc__)
		return

	wd = window

	# mean of idata
	mu = np.nanmean(idata)

	# standard deviation
	sig = np.std(idata)

	# normalized time series 
	Z = (idata[:] - mu)/sig

	# total length
	N = 2.*wd

	Tk = np.zeros(len(idata))
	for i in range(1,len(idata)):
		# first few months
		if i <= wd:
			Tk[i] = (len(Z[0:i])*(np.mean(Z[0:i]) - np.mean(Z[0:i+wd]))**2. + N/2.*(np.mean(Z[i:i+wd]) - np.mean(Z[0:i+wd]))**2.)/np.std(Z[0:i+wd])**1.
		# last few months
		elif i >= len(idata)-wd:
			Tk[i] = (N/2.*(np.mean(Z[i-wd:i]) - np.mean(Z[i-wd:]))**2. + len(Z[i:])*(np.mean(Z[i:]) - np.mean(Z[i-wd:]))**2.)/np.std(Z[i-wd:])**1.
		# everything else
		else:
			Tk[i] = (N/2.*(np.mean(Z[i-wd:i]) - np.mean(Z[i-wd:i+wd]))**2. + N/2.*(np.mean(Z[i:i+wd]) - np.mean(Z[i-wd:i+wd]))**2.)/np.std(Z[i-wd:i+wd])**1.

	return Tk
	

def Apply_MovingAverage(idata=[],window=12):
	"""
	Input:  
		idata (list of floats) input timeseries
		window (int) moving average window, default is 12 (i.e., total window size -> +-window/2)
	Output: 
		odata (list of floats) smoothed timeseries
		
	Notes:
		* The first and last window/2-timesteps are masked with NaN.
	"""

	import numpy as np

	if len(idata) == 0:
		print(f'In routine: {_currentroutine}')
		print(Apply_MovingAverage.__doc__)
		return

	rmw = np.ones(window)
	ndata = len(idata)
	odata = np.ma.zeros(ndata)

	for i in range(ndata):					# loop over all DATA 
		switch0 = 1					# switch: 
		for j in range(len(rmw)):			# loop over all WEIGHTS
			rmi = i - (len(rmw)-1)//2 + j		# running mean index: data index (i) minus half number of weights plus weight index (j)
											# asymmetrisch: dh bei len(rmw)==12 werden die 5 Monate vor und 6 Monate nach Zielmonate verw.
			if rmi < 0: continue			# if index is smaller than 0, continue to next weight. important for left boundary
			if rmi > ndata-1: continue		# important for right boundary
			odata[i] = odata[i] + idata[rmi]*rmw[j]
			index1 = j							# index of the last weight that is used
			if switch0 == 1:					
				index0 = j						# index of the first weight
				switch0 = 0						# turn off switch
		odata[i] = odata[i]/np.sum(rmw[index0:index1+1])	# divide by sum of the weights

	MC = np.ma.array([0],mask=[True])[0]
	odata[:window//2] = MC #np.NaN
	odata[-window//2:] = MC #np.NaN
	return odata


def Compute_Tendency(idata=[],DateArr=[]):
	"""
	This function computes the temporal tendency of a field (e.g., OHCT).
	
	Input:
		idata (list of floats) input timeseries.
		DateArr (list) list of datetime objects.
		
	Output:
		tendencies
		
	Notes:
		* The tendency of the first and last time step are computed with forward and backward finite differences.
		* At all other steps central finite differences are used.
	"""

	import numpy as np
	
	if len(idata) == 0:
		print(f'In routine: {_currentroutine}')
		print(Compute_Tendency.__doc__)
		return
	
	idata = np.ma.array(idata)
	
	odata = np.zeros(len(idata))
	
	# first timestep
	odata[0] = (idata[1] - idata[0])/(DateArr[1]-DateArr[0]).total_seconds()
	
	# last timestep
	odata[-1] = (idata[-1] - idata[-2])/(DateArr[-1]-DateArr[-2]).total_seconds()
		
	# all other steps
	for i in range(1,len(idata)-1):
		odata[i] = (idata[i+1] - idata[i-1])/(DateArr[i+1]-DateArr[i-1]).total_seconds()

	return odata	


def Compute_Correlation(idata1=[], idata2=[], detrend=False, return_value=False):
	"""
	input:
		idata1 (list of floats) input timeseries 1
		idata2 (list of floats) input timeseries 2
		detrend (bool) detrend input timeseries
		return_value (bool) return (or print) correlation coeffient
	"""	
	
	import numpy as np
	import scipy.signal as signal 
	
	if len(idata1) == 0 or len(idata2) == 0:
		print(f'In routine: {_currentroutine}')
		print(Compute_Correlation.__doc__)
		return
	
	if detrend:
		idata1_mod = signal.detrend(idata1[:])
		idata2_mod = signal.detrend(idata2[:])
	else:
		idata1_mod = idata1[:]
		idata2_mod = idata2[:]
	
	if return_value == False:
		return print(f'Pearson correlaion = {np.corrcoef(idata1_mod,idata2_mod)[0,1]:.3f}')
	else:
		return np.corrcoef(idata1_mod,idata2_mod)[0,1]


def Compute_AnomalyClimatology(DateArr=[],idata=[],reftime = [None]):
	"""
	Input: 
		DateArr (list of datetime objects) time 
		idata (list of floats) input timeseries
		reftime (list of 4 int) reference time in the format [start_year, start_month, end_year, end_month]
	
	Output: 
		Anomaly, Climatology
	"""

	import numpy as np
	import datetime

	if len(DateArr) == 0 or len(idata) == 0:
		print(f'In routine: {_currentroutine}')
		print(Compute_AnomalyClimatology.__doc__)
		return

	ndata = len(idata)
	odata_anom = np.ma.zeros(ndata)
	odata_clim = np.ma.zeros(12)
	iorder = np.arange(12)

	if reftime[0] == None:
		if DateArr[0].month != 1:
			iorder = np.roll(iorder,-(12-DateArr[0].month+1))

		for i in range(12):
			odata_anom[i::12] = idata[i::12] -  np.nanmean(idata[i::12])
			odata_clim[i] = np.nanmean(idata[iorder[i]::12])		
	else:
		ref_clim = np.ma.zeros(12)
		ref_sdate = datetime.datetime(reftime[0],reftime[1],DateArr[0].day)
		ref_edate = datetime.datetime(reftime[2],reftime[3],DateArr[0].day)

		si = 0
		ei = -1
		for i in range(len(idata)):
			if DateArr[i] == ref_sdate:
				si = i
			if DateArr[i] == ref_edate:
				ei = i+1

		ref_data = idata[si:ei]
		reftime = DateArr[si:ei]
				

		if reftime[0].month != 1:
			iorder = np.roll(iorder,-(12-reftime[0].month+1))	
		for i in range(12):
			ref_clim[i] = np.nanmean(ref_data[iorder[i]::12])	
		
		for i in range(len(idata)):
			odata_anom[i] = idata[i] - ref_clim[DateArr[i].month-1]

		odata_clim = np.copy(ref_clim)
		
	return odata_anom, odata_clim


def Compute_ConfidenceIntervals(idata=[], DateArr=[], CL = 0.95, ac='ar1', ano = True, print_vals = True):
	"""
	This function computes the trend and corresponding confidence intervals (CI) of an anomaly timeseries.

	Input:
		idata (list of floats) input timeseries
		DateArr (list of datetime objects) time
		CL (float) confidence level, default is 0.95
		ac (str) defines how the effective sample size is computed, possible values are 'ar1' (default) and 'full'
		ano (bool) True if idata are anomalies, default is True
	
	Output:
		Upper CI, trend, lower CI
		
	Notes:
		* This function computes confidence intervals for the trend based on anomalies using a two-tailed t-test
		* The effective sample size (auto-correlation) is considered for computing degrees of freedom and residual variance
		* T-test checks whether the trend deviates from 0 (no trend) or not
		* Default confidence level (CL) is 95 %
	"""
	
	import numpy as np
	import scipy.stats as stats
	from sklearn import linear_model
    
	if len(idata) == 0 or len(DateArr) == 0:
    		print(f'In routine: {_currentroutine}')
    		print(Compute_ConfidenceIntervals.__doc__)
    		return
    
	if ano == True:
		idata_anom = np.copy(idata)
	else:
		idata_anom, _ = Compute_AnomalyClimatology(DateArr,idata)
    
	xdata = np.arange(len(DateArr))  

	model_ols =  linear_model.LinearRegression()
	model_ols.fit(xdata.reshape(-1, 1),idata_anom) 
	idata_anom_trend = model_ols.coef_[0]   

	y_fit = model_ols.intercept_ + xdata*idata_anom_trend

	### regression residual
	regression_residual = idata_anom - y_fit

	r = np.corrcoef(regression_residual[1:],regression_residual[:-1])[0,1]
	Ntrue = len(idata)
	Neff = EffectiveSampleSize(regression_residual,ac=ac)

	### Standard error
	residual_variance = 1/(Neff)*np.sum(regression_residual**2.)
	Standard_Error = np.sqrt(residual_variance/(np.sum((xdata - np.mean(xdata))**2.)))   

	alpha = 1 - CL
	df = Neff - 2
	tvalue = stats.t.ppf(1. - alpha/2.,df)

	beta_0 = 0. # i.e. the null hypothesis is that x and y are uncorrelated -> no trend
	tstatistic = np.abs((idata_anom_trend - beta_0)/Standard_Error)

	trend_signif = ('NOT ' if tstatistic < tvalue else '')

	### ALTERNATIVE COMPUTATION
	r_time = np.corrcoef(xdata,idata_anom)[0,1]
	tstatistic_estimate = np.abs(r_time*np.sqrt(Neff)/np.sqrt(1. - r_time**2))
	CI_lower_estimate = stats.norm.ppf(alpha/2., loc=idata_anom_trend, scale=Standard_Error)
	CI_upper_estimate = stats.norm.ppf(1. - alpha/2., loc=idata_anom_trend, scale=Standard_Error)

	###########################

	CI_upper = idata_anom_trend + tvalue*Standard_Error
	CI_lower = idata_anom_trend - tvalue*Standard_Error

	if print_vals: print(f' -- (t value,t stat., t stat. est.) = ({tvalue:.1f},{tstatistic:.1f},{tstatistic_estimate:.1f})')
	if print_vals: print(f' --    Confidence intervals = ({CI_upper:.3f},{CI_lower:.3f})')
	if print_vals: print(f' -- Alternative estimate CI = ({CI_upper_estimate:.3f},{CI_lower_estimate:.3f})\n\n')        
	if print_vals: print(f'# Trend is {trend_signif}significant at CL = {CL*100:.1f} %.')        
    
	return np.array([CI_upper,idata_anom_trend,CI_lower])





