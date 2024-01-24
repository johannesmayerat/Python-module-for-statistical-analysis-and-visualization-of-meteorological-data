# Python functions for gridded data.
# 
# AUTHOR: Johannes Mayer (johannes.mayer@univie.ac.at, johannes.mayer.private@gmail.com)
# 
# VERSION: 2023-11-07
# 
# FUNCTIONS:
# 	ReadData
# 	PlotData
# 	MaskData
# 	SaveData
# 	CutToDateRange
# 	CreateDateArray
# 	Compute_Statistics
# 	Compute_PatternCorrelation
# 	Compute_Tendency
# 	Compute_MeanOverTime
# 	Compute_AnomalyClimatology
# 	Compute_TrendSignificance
# 	Apply_PredefinedMask
# 	_mila
# 	_milo
# 	_mclo
# 	_mcla
# 	_macola
# 	_get_coords
# 	_FindDateIndex
# 	_CalculateAreas
# 	_DefineMask
# 	_LineSelection
# 	

from shared import *
_currentroutine = 'gridded.py'

_Path_Mask = '/PATH/TO/LS-MASK/'
_default_opath_plot = './'
_default_opath_data = './'


def ReadData(ifile = None,varname = None):
	"""
	Input: 
	 	ifile (str) path to grib or netcdf file
	 	varname (str; optional) only for netcdf; selects the variable if netcdf-file contains multiple variables.
	Output: 
		data, DateArr, lats, lons, gridType, level, mask, units
		
	Notes:
		* GRIB data must be 2- or 3-dimensional, netcdf data can also have more than 3 dimensions.
	"""
	
	import numpy as np
	import eccodes as ecc
	import datetime
	from netCDF4 import Dataset
	from inspect import currentframe, getframeinfo
	import spharm
	from sys import exit

	if ifile == None:
		print(f'In routine: {_currentroutine}')
		print(ReadData.__doc__)
		return

##############
# GRIB FILE 
##############

	if ifile[-5:] == '.grib': 

		# DEFINE LATS LONS
		iff = open(ifile,'r')
		nmessage = ecc.eccodes.codes_count_in_file(iff)
		igrib = ecc.codes_grib_new_from_file(iff)
		nlat = ecc.codes_get(igrib,'Nj')
		nlon = ecc.codes_get(igrib,'Ni')
		latOFGP = ecc.codes_get(igrib,'latitudeOfFirstGridPointInDegrees')
		latOLGP = ecc.codes_get(igrib,'latitudeOfLastGridPointInDegrees')
		lonOFGP = ecc.codes_get(igrib,'longitudeOfFirstGridPointInDegrees')
		lonOLGP = ecc.codes_get(igrib,'longitudeOfLastGridPointInDegrees')
		GridType = ecc.codes_get(igrib,'gridType')
		param = ecc.codes_get(igrib,'paramId')
		varname = ecc.codes_get(igrib,'name')
		units = ecc.codes_get(igrib,'units')
		ecc.codes_release(igrib)
		iff.close()
		
		
		if GridType == 'reduced_gg':
			frameinfo = getframeinfo(currentframe())
			print(f'ERROR! Module_Read: Data on reduced Gaussian grid (Module_Read.py:line {frameinfo.lineno+2}).')
			exit()
		elif GridType not in ['regular_gg','regular_ll']:
			frameinfo = getframeinfo(currentframe())
			print(f'ERROR! Module_Read: Unknown grid (Module_Read.py:line {frameinfo.lineno+2}).')
			exit()

		if GridType == 'regular_gg':
			lats = spharm.spharm.gaussian_lats_wts(nlat)[0]
		else:
			GridResol = 180./nlat
			lats = [90.-GridResol/2.-GridResol*i for i in range(nlat)]
		
		if nlon > 9999: nlon = 1280

		if latOLGP > latOFGP: 
			tmp = latOFGP
			latOFGP = latOLGP
			latOLGP = tmp
		if lonOFGP > lonOLGP: lonOFGP = lonOFGP - 360.

		dlon = (lonOLGP-lonOFGP)/(nlon-1)  
		lons = np.zeros([nlon])
		for i in range(nlon):
			lons[i] = lonOFGP + i*dlon 

		if round(lats[0],2) != round(latOFGP,2) or round(lats[-1],2) != round(latOLGP,2):
			print(f'{round(lats[0],2)} != {round(latOFGP,2)} ?') 
			print(f'{round(lats[-1],2)} != {round(latOLGP,2)} ?') 
			frameinfo = getframeinfo(currentframe())
			print(f'ERROR! Module_Read: Wrong grid (Module_Read.py:line {frameinfo.lineno+2}).')
			exit()

		if round(lons[-1],2) != round(lonOLGP,2):
			print(f'{round(lons[-1],2)} != {round(lonOLGP,2)}') 
			frameinfo = getframeinfo(currentframe())
			print(f'ERROR! Module_Read: Wrong grid (Module_Read.py:line {frameinfo.lineno+2}).')
			exit()

		#lons = lons[:] - 180.	


		# READ IN DATA
		dataset = np.zeros([nlat,nlon,nmessage])
		mask = np.zeros([nlat,nlon,nmessage])
		time = []		
		iff = open(ifile)
		for i in range(nmessage):
			igrib = ecc.codes_grib_new_from_file(iff)

			if igrib is None:
				break

			dataset[:,:,i] = ecc.codes_get_values(igrib).reshape([nlat,nlon])

			time += [datetime.datetime.strptime(str(ecc.codes_get(igrib,'dataDate')),'%Y%m%d')]
			ecc.codes_release(igrib)
		iff.close()

		if all(i == time[0] for i in time):
			level = np.arange(nmessage)
			time = [time[0]]
			dim3 = 'level'
		else:
			level = [0] #np.zeros(nmessage)
			dim3 = 'time'

		if nmessage == 1:
			dataset = dataset[:,:,0]
			mask = mask[:,:,0]

		dataset = np.ma.array(dataset)
		
################
# NETCDF FILE
################

	elif ifile[-3:] == '.nc':
		#print('NC')

		idata = Dataset(ifile, mode='r')

		List_LatsName = ['lat','latitude','ncl2'] # insert all possible lat dimension names
		List_LonsName = ['lon','longitude','ncl3']
		List_TimeName = ['time','Time','time_counter','ncl1']
		List_LevelName = ['level','lev']
		lats = []		
		lons = []
		time_unf = []
		level = [0]
		for i in List_LatsName:
			try:
				lats = idata.variables[i][:]
				latsVarName = i
				break
			except:
				latsVarName = None

		for i in List_LonsName:
			try:
				lons = idata.variables[i][:]
				lonsVarName = i
				break
			except:
				lonsVarName = None

		for i in List_TimeName:
			try:
				time_unf = idata.variables[i][:]
				timeVarName = i
				break
			except:
				timeVarName = None

		for i in List_LevelName:
			try:
				level = idata.variables[i][:]
				levelVarName = i
				break
			except:
				levelVarName = None	
	
		if lats == [] or lons == []:
			frameinfo = getframeinfo(currentframe())
			print(f'ERROR! Input dimensions not found (Module_Read.py:line {frameinfo.lineno+2})')
			print(idata.variables)
			exit()		

		if all([abs(lats[i] - lats[i+1]) == abs(lats[1] - lats[0]) for i in range(len(lats)-1)]):
			GridType = 'regular_ll'
		else:
			GridType = 'regular_gg'

		time = []
		if time_unf == []:
			print('No time information given. Select format and time.')
			format_time = input("Format (e.g., '%Y-%m-%d'): ")
			time += [datetime.datetime.strptime(input('Time: '), format_time)]
		else:
			dataunits = idata.variables[timeVarName].units.split(' ')

			if 'seconds' in dataunits:
				refdate = datetime.datetime.strptime(dataunits[2].replace('T00:00:00+00:00',''),'%Y-%m-%d')
				time += [refdate + datetime.timedelta(seconds=int(i)) for i in time_unf]	
			elif 'hours' in dataunits:
				refdate = datetime.datetime.strptime(dataunits[2],'%Y-%m-%d')
				time += [refdate + datetime.timedelta(hours=int(i)) for i in time_unf]	
			elif 'days' in dataunits:
				refdate = datetime.datetime.strptime(dataunits[2],'%Y-%m-%d')
				time += [refdate + datetime.timedelta(days=int(i)) for i in time_unf]	
			elif 'YYYYMM' in dataunits:
				time += [datetime.datetime.strptime(str(int(i)).replace('.',''),'%Y%m') for i in time_unf]
			else:
				frameinfo = getframeinfo(currentframe())
				print(f'ERROR! Unknown unit format (Module_Read.py:line {frameinfo.lineno+2}).')
				exit()

		List_Vars = idata.variables.keys()
		List_Dims = idata.dimensions.keys()
		if len(List_Vars) == len(List_Dims)+1:
			for i in List_Vars:
				if i not in List_Dims: varname = i
		else:
			if varname == None:
				for i in List_Vars: print(' ',i)
				varname = input('Multiple variables available. Select: ')
	
		while varname not in List_Vars:
				for i in List_Vars: print(' ',i)
				varname = input(f'{varname} not found in input. Select: ')

		dataset = np.squeeze(idata.variables[varname][:])

		dataset.mask = [False]
		dataset.mask[dataset > 1e30] = [True]
		dataset.mask[np.isnan(dataset)] = [True]

		datadims = []
		for i in range(len(idata.variables[varname].get_dims())):
			if idata.variables[varname].get_dims()[i].size != 1:
				datadims += [idata.variables[varname].get_dims()[i]]


		while len(dataset.shape) >= 4:
			print(f'Dataset "{varname}" has more than 3 dimensions:')
			cnt_dims = 0
			#datadims = []
			allowed_dims = []
			for i in range(len(datadims)):
				if datadims[i].size != 1:
					if datadims[i].name in [latsVarName,lonsVarName]:
						print(f' (X) {datadims[i].name} (size = {datadims[i].size})')
					else:
						print(f' ({cnt_dims}) {datadims[i].name} (size = {datadims[i].size})')
						allowed_dims += [cnt_dims]
					cnt_dims += 1
			print('Select dimension and index.')
			sel_dim = int(input('Select dimension: '))
			while sel_dim not in allowed_dims:
				print(f'Not possible to select dimension {sel_dim}.')
				sel_dim = int(input('Dimension: '))

			for i in range(datadims[sel_dim].size):
				print(f' ({i}) {datadims[sel_dim].name} = {idata.variables[datadims[sel_dim].name][i]}')

			print(f'Index range: 0 .. {datadims[sel_dim].size-1}')
			sel_ind = int(input('Select index: '))
			while sel_ind < 0 or sel_ind > datadims[sel_dim].size-1:
				print(f'Not possible to select index {sel_ind}.')
				sel_ind = int(input('Index: '))

			if datadims[sel_dim].name == timeVarName:
				time = [time[sel_ind]]
			if datadims[sel_dim].name == levelVarName:
				level = [level[sel_ind]]

			datadims.pop(sel_dim)
			dataset = dataset.ma.take(sel_ind,axis=sel_dim)
			
		if len(dataset.shape) == 3:
			ddNames = []
			for i in range(len(datadims)):
				ddNames += [datadims[i].name]
				if datadims[i].name not in [latsVarName,lonsVarName]: 
					dim3i = i
					dim3 = datadims[i].name

			#print(ddNames)

			if ddNames[0] != latsVarName or ddNames[1] != lonsVarName:
				Permutation = (ddNames.index(latsVarName),ddNames.index(lonsVarName),ddNames.index(ddNames[dim3i]))
				print(Permutation)
				dataset = np.ma.transpose(dataset, Permutation)
			del ddNames
		else:
			dim3 = None

		if np.ma.is_masked(dataset):
			mask = dataset.mask
		else:
			mask = np.zeros_like(dataset)
			dataset = np.ma.array(dataset)
		
		#units = None
		List_UnitName = ['units','Units','unit','Units']
		for i in List_UnitName:
			try:
				units = getattr(idata.variables[varname],i)
				break
			except:
				pass

		if units == None:
			print('No unit argument given. Using [W m**-2]')
			units = 'W m**-2'
		else:
			print(f'Input units found: {units}')
			
		param = 9999

	else:
		frameinfo = getframeinfo(currentframe())
		print(f'ERROR! Unknown file format (Module_Read.py:line {frameinfo.lineno+2}).')
		exit()

	if lats[0] < 0.: 
		lats = np.flip(lats)
		dataset = np.flip(dataset,axis=0)
		mask = np.flip(mask,axis=0)

	if dataset.size != dataset.mask.size:
		if dataset.mask.size == 1:
			dataset.mask = [dataset.mask]
		else:
			frameinfo = getframeinfo(currentframe())
			print(f'ERROR! Unsuitable data mask (Module_Read.py:line {frameinfo.lineno+2}).')
			exit()

	return dataset, time, lats, lons, GridType, level, mask, units



def PlotData(idata=[],lats=[],lons=[],\
		units=None, pltrange=None, zmmax=None, nozm=False, nocbar=False, cmap='default',\
		plot_shading=None, plot_subregion=None, plot_vector=None, plot_rect=None, plot_contour=None,\
		text_corner=None, text_left=None, text_title=None,\
		labelleft=True, labelbottom=True, labeltop=False, labelright=False,\
		opath=_default_opath_plot, savefig=None): # camp='RdBu_r'
		
	"""
	Necessary arguments:
		idata	(2D array) input data with dimensions (nlats,nlons)
		lats	(list of floats) latitudes
		lons	(list of floats) longitudes
		units	(str), '**' will be translated to exponent
		
	Optional arguments:
		pltrange (list of floats) min and max values of the plot 
		zmmax 	 (float) zonal mean maximum (symmetric)
		nozm	 (bool) hide zonal mean plot
		nocbar	 (bool) hide color bar
		cmap	 (matplotlib.colors.LinearSegmentedColormap object) color map
		
		plot_shading   (2D array) plot shading array, dimenions of idata
		plot_subregion (list of str) plot sub-region (e.g., ['00N', '90N', '90W', '30E'])
		plot_vector    (list) draws vectors over the current plot, input is a list with meridional and zonal winds and scale
		plot_rect      (list of lists of str) draws rectangles over given area, e.g., [['76N','68N','2E','10E'],['20N','30N','50W','40W']]
		plot_contour   (list) draws contour lines at given values
		
		text_corner (str) text in the left corner inside the plot
		text_left   (str) text in the left corner outside the plot
		text_title  (str) title text of the plot; if None, then MEAN and RMS are shown.
		
		labelleft, labelbottom, labeltop, labelright (bool) turn on/off labels
		opath   (str) output path
		savefig (str) saves plot with 'savefig' as filename (e.g., 'TEDIV.png')

	"""
	
	import matplotlib
	import numpy as np
	import matplotlib.pylab as plt
	import matplotlib.patches as mpatches
	from matplotlib.ticker import MultipleLocator
	from matplotlib.ticker import FormatStrFormatter
	from matplotlib.colors import LinearSegmentedColormap	
	import cartopy.crs as ccrs
	from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter	
	from sys import exit
	
	# ignore deprecation warning:  "ShapelyDeprecationWarning: __len__ for multi-part geometries is deprecated and will be removed in Shapely 2.0. 
	# Check the length of the `geoms` property instead to get the  number of parts of a multi-part geometry."
	import warnings
	from shapely.errors import ShapelyDeprecationWarning
	warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning) 

	from shared import _rcparams
	fs_all, fw_all = _rcparams()

	if len(idata) == 0 or len(lats) == 0 or len(lons) == 0 or units == None:
		print(f'In routine: {_currentroutine}')
		print(PlotData.__doc__)
		return

	if cmap=='default': cmap = matplotlib.colors.LinearSegmentedColormap.from_list('custom',['#1f77b4','white','firebrick']) # '#d62728'


	##  PREPARE TITLE
	if text_title != None and len(text_title) > 50:

		# search spaces in text_title string
		ispace = []
		for i in range(len(text_title)):
			if text_title[i] == ' ': ispace += [i]

		# compute central index 
		imid = len(text_title)//2

		# evaluate space closest to the center of the string
		spi = ispace[np.argmin(np.abs(np.array(ispace)-imid))]

		text_title = text_title[:spi]+'\n'+text_title[spi:]
	elif text_title == None:
		data_mean,data_rms, *_ = Compute_Statistics(idata,lats,lons)
		if all([ i < 1e-2 for i in [data_mean,data_rms]]):
			text_title = f'MEAN = {data_mean:.2e}, RMS = {data_rms:.2e}'
		else:
			text_title = f'MEAN = {data_mean:.2f}, RMS = {data_rms:.2f}'
			
	if pltrange == None: pltrange = [-100,100]

	## PREPARE UNITS
	PlotUnits = ''
	for i in units.split(' '): # 
		end_str = ' '
		if '**' in i: end_str = '}$ '		
		PlotUnits += i.replace('**',r'$^{') + end_str

	DG = r'°' #  r'$^{\circ}$'
	DGW = r'°W'
	DGE = r'°E' # r'$^{\circ}$E'
	DGS = r'°S'
	DGN = r'°N'
	XLBL = np.array(['180'+DG,'120'+DGW,'60'+DGW,'0'+DG,'60'+DGE,'120'+DGE,'180'+DG])
	XTIC = np.array([-179.99999, -120, -60, 0, 60, 120, 179.99999])
	YLBL = np.array(['90'+DGS,'60'+DGS,'30'+DGS,'0'+DG,'30'+DGN,'60'+DGN,'90'+DGN])
	YTIC = np.array([-90, -60, -30, 0, 30, 60, 90])


	fig = plt.figure(figsize=(19,10))
	fig.subplots_adjust(hspace=0, wspace=0.0)

	if nozm:
		sel_colspan = 10
	else:
		sel_colspan = 9

	ax = plt.subplot2grid((5,10),(0,0), colspan=sel_colspan, rowspan=5, projection=ccrs.PlateCarree(central_longitude=0))
	ax.coastlines()
	ax.background_patch.set_facecolor('grey') # 'grey' 'silver'


	## PLOT SUB REGION
	if plot_subregion == None:
		ax.set_global()
	else:		
		lat0_plt,lat1_plt,lon0,lon1 = _get_coords(plot_subregion)	
		if lon0 > 180 : lon0 -= 360
		if lon1 > 180 : lon1 -= 360
		
		if lon1 < lon0:
			plt.xlim(lon1,lon0)
		else:
			plt.xlim(lon0,lon1)
		plt.ylim(lat1_plt,lat0_plt)
				
	if plot_subregion == None:
		ax.set_xticks(XTIC, crs=ccrs.PlateCarree())
		ax.set_xticklabels(XLBL)
		ax.set_yticks(YTIC, crs=ccrs.PlateCarree())
		ax.yaxis.set_major_formatter(LatitudeFormatter())
		ax.xaxis.set_minor_locator(MultipleLocator(15))
		ax.yaxis.set_minor_locator(MultipleLocator(10))
	else:
		ylim = ax.get_ylim()
		xlim = ax.get_xlim()
		dy = np.abs(ylim[1]-ylim[0])/3
		dx = np.abs(xlim[1]-xlim[0])/3
		x_arr = np.arange(xlim[0],xlim[1]+dx,dx)
		y_arr = np.arange(ylim[0],ylim[1]+dy,dy)
		ax.set_xticks(x_arr)
		ax.set_yticks(y_arr)
		
		ax.yaxis.set_minor_locator(MultipleLocator(10))
		
		x_arr_lbl = []
		y_arr_lbl = []
		for i in range(len(x_arr)):
			x_arr_lbl += [str(int(abs(x_arr[i])))+(DGW if x_arr[i] < 0 else DGE if x_arr[i] > 0 else DG)]
		for i in range(len(y_arr)):
			y_arr_lbl += [str(int(abs(y_arr[i])))+(DGS if y_arr[i] < 0 else DGN if y_arr[i] > 0 else DG)]				
			
		ax.set_xticklabels(x_arr_lbl)
		ax.set_yticklabels(y_arr_lbl)			
		
		
	ax.tick_params(axis='both', which='major', pad=6, labelsize=fs_all,direction='out')
	ax.tick_params(axis='both', which='minor', length=6, direction='out')
	
	
	## PLOT DEFAULT
	cs = ax.pcolormesh(lons, lats, idata, cmap=cmap, vmin=pltrange[0], vmax=pltrange[1], linewidth=0, rasterized=True, shading='auto')
	
	
	## PLOT VECTOR FIELD 
	if plot_vector != None:

		udata = plot_vector[0]
		vdata = plot_vector[1]
	
		ds = 10
		if len(plot_vector) > 2:
			if plot_vector[2] == None:
				quiver_scale = ds//2
			else:
				quiver_scale = plot_vector[2]
		else:
			quiver_scale = ds//2
			
		#print(f'u max/min: {udata.max():.4f}/{udata.min():.4f}')
		#print(f'v max/min: {vdata.max():.4f}/{vdata.min():.4f}')
		
		plt.quiver(lons[::ds],lats[::ds],udata[::ds,::ds],vdata[::ds,::ds],scale=quiver_scale,width=0.001)
		
	
	## PLOT CONTOUR LINES
	if plot_contour != None:
		plt.contour(lons, lats, idata, levels=plot_contour,colors=['k'],linestyles=[':'],linewidths=1)
		
	
	## PLOT SHADING
	if not isinstance(plot_shading,type(None)):
		if not isinstance(plot_shading,np.ndarray):
			print('ERROR! Unknown data type for plot_shading.')
			exit()
		
		plt.contourf(lons, lats, plot_shading, 1, hatches=['','x'],extend='lower',alpha=0) # hatches=['','..']
	
	## PLOT RECTANGLES
	if plot_rect != None:
		for i in range(len(plot_rect)):
			lat0,lat1,lon0,lon1 = _get_coords(plot_rect[i])
			ax.add_patch(mpatches.Rectangle(xy=[lon0,lat1], height=np.abs(lat0-lat1), width=np.abs(lon1-lon0), edgecolor='k', facecolor='none', ls='-', lw=3, alpha=1.0, transform=ccrs.PlateCarree()))

	## TITLE TEXT
	plt.title(text_title,loc='center', fontsize=fs_all, weight=fw_all, pad=10)	
	
	## CORNER TEXT (inside the plot)
	if text_corner != None:
		ax.text(0.02,0.97,text_corner,fontsize=fs_all, bbox=dict(facecolor='white',edgecolor='black'),transform=ax.transAxes, ha='left', va='top')
	
	## LEFT TITLE TEXT
	if text_left != None:
		plt.title(text_left, loc='left', fontsize=fs_all, weight=fw_all, pad=10 )
		
	
	## COLOR BAR
	if not nocbar: 
		if abs(pltrange[1]) < 1e-1:
			cbar_format = '%.1e'
		elif abs(pltrange[1]) <= 1e-0:
			cbar_format = '%.1f'
		elif abs(pltrange[1]) >= 1e5:
			cbar_format = '%.1e'
		else: 
			cbar_format = '%.0f'

		cbar = plt.colorbar(cs,ax=ax,orientation='horizontal', pad=0.09, shrink=0.55,aspect=20, anchor=(0.5,0.0),format=cbar_format) 
		cbar.set_label(PlotUnits,fontsize=fs_all,weight=fw_all)
		cbar.ax.tick_params(labelsize=fs_all)
		cbar.ax.locator_params(nbins=5)
		
		if plot_contour != None:
			contour_ext = [[i,i] for i in plot_contour]
			for i in range(len(plot_contour)): 
				cbar.ax.plot(contour_ext[i],([[-1000, 1000]]*len(plot_contour))[i],'k:')
				
				
	## ZONAL MEAN
	if nozm:				 
		plt.tick_params(axis='y', which='both', labelleft=labelleft, labelright=labelright)
		plt.tick_params(axis='x', which='both', labeltop=labeltop, labelbottom=labelbottom)
		
	else:
		## horizonale posi, vertikale posi, breite, höhe
		ax2 = plt.Axes(fig, [ax.get_position().x1-0., ax.get_position().y0, .07, ax.get_position().y1-ax.get_position().y0])	
		fig.add_axes(ax2)

		## Compute zonal mean
		ZonalMean = np.flip(np.nanmean(idata,axis=1),axis=0)

		## Define y grid
		ygrid = np.linspace(lats[-1],lats[0],len(lats))

		## Compute limits for plot
		xlimit = np.nanmax(abs(ZonalMean[5:len(lats)-5]))

		## Preamble
		plt.tick_params(labelleft=False, labelright=True, labeltop=True, labelbottom=False)
		plt.yticks(YTIC,YLBL,rotation=0,fontsize=fs_all,weight=fw_all)
		plt.xticks(rotation=90,fontsize=fs_all-5,weight=fw_all)

		if zmmax == None:
			plt.xlim(-xlimit*1.1,xlimit*1.1)
		else:
			plt.xlim(-zmmax,zmmax)
			
		if plot_subregion == None:
			plt.ylim(-90,90)
		else:
			plt.ylim(lat1_plt,lat0_plt)
			ax2.set_yticks(y_arr)
			ax2.set_yticklabels(y_arr_lbl)
			
			
		plt.grid(linestyle='--')
		ax2.yaxis.set_minor_locator(MultipleLocator(10))
		if abs(pltrange[1]) < 1e-1:
			ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1e'))
			plt.xticks(fontsize=fs_all-10,rotation=90.)

		## ZONAL MEAN PLOT
		ax2.plot(ZonalMean,ygrid,"k-")
		ax2.fill_betweenx(ygrid,ZonalMean,np.zeros(len(ygrid)), where=ZonalMean <= 0., color='lightblue',alpha=0.6)
		ax2.fill_betweenx(ygrid,ZonalMean,np.zeros(len(ygrid)), where=ZonalMean >= 0., color='lightcoral',alpha=0.6)

	if savefig != None:
		if '.p' not in savefig: savefig += '.png'
		
		
		plt.savefig(opath+savefig,dpi=300,bbox_inches='tight',transparent=False,facecolor='white')
		print(f' -- Figure saved: {opath+savefig}')
	return


def MaskData(idata=[], ilats=[], ilons=[],\
		amask=False, vmask=False, lmask=False, pmask=False,\
		gridType='regular_gg',print_info=False):
	"""
	This function applies certain types of masks on the data.
	
	Input:
		idata (2D/3D array) input data
		ilats (list of floats) latitudes in degree
		ilons (list of floats) longitudes in degree
		amask (str) area mask, e.g., '00N 90N 90W 30E'
		vmask (str) value mask, first two characters must be alphabetic symbols, e.g., 'gt1e10'
		lmask (str) latitude mask, e.g., '60E', '70N', '40N60S' 
		pmask (str) predefined masks, e.g., 'land', 'ocean',... (see predefined masks in function _DefineMask)
		gridType (str) grid type of the data, can be 'regular_gg' or 'regular_ll'
		
	Output:
		masked_data (same as idata but masked)
		
	Notes:
		lmask:  '60E' is a symmetric mask around the equator where everything north and south of 60° is masked.
			'60N' masks everything north of 60°N.
			'40N60S' masks everything north of 40°N and south of 60°S.
			
		vmask:  Possible arguments are 'ge', 'le', 'gt', 'lt' (greater/less equal/than) + a numeric value, together as string, e.g., 'gt1e10'.
			If vmask is just a numeric value (which must still be parsed as string), every grid point with exactly this value is masked.
		
		pmask:	'Indic', 'Atlantic', 'Pacific', 'AtlanticArctic', 'ARCGATE2FRAM', 'ARCGATE', 'Atlantic2Fram',
			'DomainN+S', 'DomainN', 'DomainS', 'SAMBA2GSR', 'SAMBA2Bering', 'SAMBA2RAPID', 'Gulf', 
			'Africa', 'Antarctica', 'SouthAmerica', 'NorthAmerica', 'Australia', 'Asia', 'Europe'
	
	"""
	
	from inspect import currentframe, getframeinfo
	import numpy as np
	import re
	from sys import exit
	
	if len(idata) == 0:
		print(f'In routine: {_currentroutine}')
		print(MaskData.__doc__)
		return

	if len(idata.shape) == 2:
		idata = idata[:,:,np.newaxis]

	nlat = len(ilats)
	nlon = len(ilons)
	mask1 = np.zeros([nlat,nlon,1])

##############
# area mask	
##############

	if amask != False:
		if isinstance(amask, str):
			amask = amask.split(' ')

		if len(amask) != 4:
			frameinfo = getframeinfo(currentframe())
			print(f'ERROR!  Wrong length of amask argument ({_currentroutine}:line {frameinfo.lineno+2}).')
			exit()	
		
		lat0 = lat1 = lon0 = lon1 = None
		lat_cnt = lon_cnt = 0
		for i in range(4):
			#print(amask[i], amask[i][-1])
			if amask[i][-1] in ['N','S']:
				lat_cnt += 1
				if lat_cnt > 2:
					frameinfo = getframeinfo(currentframe())
					print(f'ERROR!  Wrong coordinates in amask argument ({_currentroutine}:line {frameinfo.lineno+2}).')
					exit()

				if lat0 == None:
					lat0 = float(amask[i][:-1])
					if amask[i][-1] == 'S': lat0 *= -1.
					continue
				if lat1 == None:
					lat1 = float(amask[i][:-1])
					if amask[i][-1] == 'S': lat1 *= -1.
					continue
			if amask[i][-1] in ['E','W']:
				lon_cnt += 1
				if lon_cnt > 2:
					frameinfo = getframeinfo(currentframe())
					print(f'ERROR!  Wrong coordinates in amask argument ({_currentroutine}:line {frameinfo.lineno+2}).')
					exit()
				if lon0 == None:
					lon0 = float(amask[i][:-1])
					if amask[i][-1] == 'W': lon0 = 360.-lon0
					continue
				if lon1 == None: 
					lon1 = float(amask[i][:-1])
					if amask[i][-1] == 'W': lon1 = 360.-lon1
					continue
	
		if lat0 == lat1 or lon0 == lon1:
			frameinfo = getframeinfo(currentframe())
			print(f'ERROR!  Wrong coordinates in amask argument ({_currentroutine}:line {frameinfo.lineno+2}).')
			exit()
		if lat0 < lat1:
			tmp = lat0
			lat0 = lat1
			lat1 = tmp
		if lon1 < lon0:
			tmp = lon0
			lon0 = lon1
			lon1 = tmp	

		#if print_info == True:
		#	print(f' -- lats: {lat0} to {lat1}')
		#	print(f' -- lons: {lon0} to {lon1}')
		
		# first mask everything
		mask1[:,:,:] = 1
		# then unmask selected region
		if lon0 <= 90 and lon1 > 240:
			# if true unmask everything outside [lon0,lon1] -> this allows to choose a box including the prime meridian
			mask1[_mcla(nlat,lat0):_mcla(nlat,lat1),:_mclo(nlon,lon0),:]  = 0
			mask1[_mcla(nlat,lat0):_mcla(nlat,lat1),_mclo(nlon,lon1):,:]  = 0

		else:
			# otherwise unmask everything inside [lon0,lon1]
			mask1[_mcla(nlat,lat0):_mcla(nlat,lat1),_mclo(nlon,lon0):_mclo(nlon,lon1),:]  = 0	
	
###############
# value mask	
###############

	if vmask != False:
		if vmask[0:2] == 'ge':
			mask1[idata[...,0] >= float(vmask[2:])] = 1
		elif vmask[0:2] == 'le':
			mask1[idata[...,0] <= float(vmask[2:])] = 1
		elif vmask[0:2] == 'gt':
			mask1[idata[...,0] > float(vmask[2:])] = 1
		elif vmask[0:2] == 'lt':
			mask1[idata[...,0] < float(vmask[2:])] = 1
		else:
			mask1[idata[...,0] == float(vmask)] = 1
			
###############
# latitude mask	
###############
	
	if lmask != False:
		# possible strings:
		# 60N, 60S, ...	: alles außer 60N bis Pol, 60S bis Pol, ...
		# 20E, 30E, ...	: alles außer 20 N/S um den Äquator, 30 N/S um den Äq., ...
		# 60N75S, ...	: asymmetrische Maske, d.h. alles zwischen 60N und 75S

		split_arr = re.split(r'([N,S,A])', lmask)[:-1]

		if not lmask[0].isnumeric():
			print(f'ERROR! Wrong invalid lmask string.')
			print(f'Examples: 60.1N, 20.4A, 60N30S')
			return

		if len(split_arr) == 2:
			lmask_val = float(split_arr[0])

			if split_arr[1] == 'N':
				mask1[_mcla(nlat,lmask_val):,:,:] = 1.
			elif split_arr[1] == 'S':
				mask1[:_mcla(nlat,-lmask_val),:,:] = 1.
			elif split_arr[1] == 'E':
				mask1[:_macola(lmask_val,ilats),:,:] = 1.		
				mask1[_macola(-lmask_val,ilats):,:,:] = 1.

		elif len(split_arr) == 4:
			lmask_val1 = float(split_arr[0])
			lmask_val2 = float(split_arr[2])

			if split_arr[1] == 'S': lmask_val1 = -lmask_val1
			if split_arr[3]	== 'S': lmask_val2 = -lmask_val2
			mask1[:_mcla(nlat,lmask_val1) ,:,:] = 1.		
			mask1[ _mcla(nlat,lmask_val2):,:,:] = 1.

		else:
			frameinfo = getframeinfo(currentframe())
			print(f'ERROR!  Wrong length of lmask argument ({_currentroutine}:line {frameinfo.lineno+2}).')
			exit()		

####################
# predefined mask	
####################

	if pmask != False: #in ['land','sea','ocean']:
		dummy_data = np.zeros([nlat,nlon])
		_, returned_mask = Apply_PredefinedMask(dummy_data,ilats,ilons,iregion=pmask,gridType=gridType,print_info=print_info)

		if np.sum(mask1[:,:,0]) == 0:
			mask1[:,:,0] = returned_mask[:,:]
		else:
			mask1[:,:,0] = np.int32(mask1[:,:,0]) | np.int32(returned_mask[:,:])	
	
	mask3D = np.repeat(mask1,idata.shape[2],axis=2)
	masked_data = np.ma.copy(idata)
	masked_data.mask[:,:,:] = mask3D[:,:,:]

	return masked_data
	
	
def SaveData(ofile=[],data=[],lats=[],lons=[],time=[],\
	varname='DATA', vn_long='DATA', units='W m**-2', reftime=None, olat_units='degrees_north', olon_units='degrees_east',\
	opath=_default_opath_data):
	"""
	This function saves the data to a NETCDF file (opath+ofile).

	Input:
		ofile (str)       output file name (path is diven by opath).
		data  (2D or 3D)  data to be saved.
		lats  (list)      latitudes in degree
		lons  (list)      longitudes in degree
		time  (list)      list of datetime objects.
		varname (str)     (optional) variable name, default is 'DATA'.
		vn_long (str)     (optional) long variable name, default is 'DATA'.
		units   (str)     (optional) units, default is 'W m**-2'.
		reftime (datetime object) (optional) reference time given as datetime object, default is 1 January 1900.
		opath   (str)     (optional) output path, default is given by _default_opath_data (see header).
	
	
	Note: for multiple variables in one file, use a list of numpy arrays, all variables must have the same dimension. 
	"""
	
	import datetime
	from netCDF4 import Dataset
	import numpy as np

	if len(ofile) == 0:
		print(f'In routine: {_currentroutine}')
		print(SaveData.__doc__)
		return
	
	if not isinstance(varname,list): varname = [varname] 
	if not isinstance(vn_long,list): vn_long = [vn_long] 
	if not isinstance(units,list): units = [units] 

	if reftime == None: reftime = datetime.datetime(1900, 1, 1, 0, 0)
	
	if isinstance(data,list):
		iter_data = len(data)
		is_3d = True if data[0].ndim == 3 else False
	elif isinstance(data,np.ndarray):
		iter_data = 1
		data = [data]
		is_3d = True if data.ndim == 3 else False
	else:
		return print('Wrong type of data array.')
	
	if is_3d and len(time) == 0:
		return print('No time coordinates given for 3D data.')
	if len(varname) != iter_data:
		return print('Not enough variable names provided.')
	if len(vn_long) != iter_data:
		return print('Not enough long variable names provided.')
	if len(units) != iter_data:
		print('Not enough units provided.')
		return 
			
	if '.nc' not in ofile: ofile += '.nc'

	f = Dataset(opath+ofile,mode='w',format='NETCDF4') 

	f.createDimension(dimname='latitude',size=len(lats))
	f.createDimension(dimname='longitude',size=len(lons))	
	if is_3d: f.createDimension(dimname='time',size=len(time))	

	olat = f.createVariable('latitude', 'f4', 'latitude')  
	olon = f.createVariable('longitude', 'f4', 'longitude')
	if is_3d: otime = f.createVariable('time', 'i4', 'time') # 'i4'

	if is_3d: 
		odata_tuple = ('latitude','longitude','time')
	else:
		odata_tuple = ('latitude','longitude')
	
	# LATITUDE COORDINATE	
	olat[:] =  np.copy(lats)
	olat.long_name = 'latitude'
	olat.units = olat_units	
	
	# LONGITUDE COORDINATE
	olon[:] =  np.copy(lons)
	olon.long_name = 'longitude'
	olon.units = olon_units	
	
	# TIME COORDINATE
	if is_3d:
		for i in range(len(time)):
			otime[i] = (time[i] - reftime).days*24 + (time[i] - reftime).seconds//3600  # time as hours since reftime
		otime.long_name = 'time'
		otime.units = 'hours since ' + datetime.datetime.strftime(reftime,'%Y-%m-%d %H:%M:%S.0')
		otime.calendar = 'gregorian'			
		
	# DATA
	odata = []
	for i in range(iter_data):
		odata.append(f.createVariable(varname[i],'f4',odata_tuple))
		odata[i][...] = data[i][...]
		odata[i].units = units[i]
		odata[i].long_name = vn_long[i]
	
	f.history =  "Created with mainr on " + datetime.datetime.today().strftime("%Y-%m-%d")
	f.close()
	print(f' -- Output written: {opath+ofile}')

	return
		

def CutToDateRange(idata = [],DateArr=[],sdate=[], edate=[], dfmt='%Y-%m-%d'):
	# idata dims: (lats, lons, time)
	"""
	This function cuts the input data and date array to the period sdate-edate (end date is included).
	
	INPUT: 
	 	idata   (3D array, with last dim=2 as time) input data
	 	DateArr (list of datetime objects) date array
	 	sdate   (str) start date in the format of dfmt
	 	edate   (str) end date in the format of dfmt
	 	dfmt    (str) datetime format, default is '%Y-%m-%d'

	OUTPUT: 
		idata:   as the input data but cut to the period sdate-edate
		DateArr: as the input date array but cut to the period sdate-edate
	"""
	
	import datetime
	import numpy as np
	
	if len(idata) == 0:
		print(f'In routine: {_currentroutine}')
		print(CutToDateRange.__doc__)
		return

	if type(sdate) == str:
		sdate = datetime.datetime.strptime(sdate,dfmt)
	if type(edate) == str:
		edate = datetime.datetime.strptime(edate,dfmt)	

	sindex = DateArr.index(sdate)
	eindex = DateArr.index(edate) + 1

	if np.array(idata).ndim == 3:
		return idata[:,:,sindex:eindex], DateArr[sindex:eindex]
	elif np.array(idata).ndim == 2:
		return idata[:,sindex:eindex], DateArr[sindex:eindex]
	else:
		return idata[sindex:eindex], DateArr[sindex:eindex]


def CreateDateArray(sdate=[],edate=[],resol='monthly', dfmt='%Y-%m-%d'):
	"""
	This function creates a list of datetime objects for the period sdate-edate (end date included).
	
	INPUT: 
		sdate (str) start date given in the format dfmt
		edate (str) end date given in the format dfmt
		resol (str) temporal resolution, possible arguments are 'hourly', 'daily', 'monthly', 'yearly'; default is 'monthly'
		dfmt  (str) datetime format, default is '%Y-%m-%d'
	OUTPUT: 
		DateArr: list of datetime objects with the resolution given by 'resol'.
	"""
	
	import datetime
	from inspect import currentframe, getframeinfo

	if len(sdate) == 0:
		print(f'In routine: {_currentroutine}')
		print(CreateDateArray.__doc__)
		return


	if type(sdate) == str:
		sdate = datetime.datetime.strptime(sdate,dfmt)
	if type(edate) == str:
		edate = datetime.datetime.strptime(edate,dfmt)	

	syear = sdate.year
	smonth = sdate.month
	sday = sdate.day
	eyear = edate.year
	emonth = edate.month
	eday = edate.day

	DateArr = []
	if resol == 'hourly':
		dt = datetime.timedelta(hours=1)
		date = sdate
		while date <= edate:
			DateArr += [date]
			date = date + dt

	elif resol == 'daily':
		for year in range(syear,eyear+1):
			DaysPerMonths = [31,28,31,30,31,30,31,31,30,31,30,31]
			if year%4 == 0: DaysPerMonths[1] = 29

			smonth_iter = 1
			emonth_iter = 13
			if year == syear:
				smonth_iter = smonth
			if year == eyear:
				emonth_iter = emonth + 1

			for month in range(smonth_iter,emonth_iter):
				sday_iter = 1
				eday_iter = DaysPerMonths[month-1] + 1
				if year == syear and month == smonth:
					sday_iter = sday
				if year == eyear and month == emonth:
					eday_iter = eday + 1

				for day in range(sday_iter,eday_iter):
					DateArr += [datetime.datetime(year,month,day)]

	elif resol == 'monthly':
		for year in range(syear,eyear+1):
			smonth_iter = 1
			emonth_iter = 13
			if year == syear:
				smonth_iter = smonth
			if year == eyear: 
				emonth_iter = emonth + 1
			for month in range(smonth_iter,emonth_iter):
				DateArr += [datetime.datetime(year,month,sday)]	
	elif resol == 'yearly':
		for year in range(syear,eyear+1):
			DateArr += [datetime.datetime(year,smonth,sday)]
	else:
		frameinfo = getframeinfo(currentframe())
		print(f'ERROR! Wrong resolution selected. Use one of the following: hourly/daily/monthly/yearly ({_currentroutine}:line {frameinfo.lineno+2}).')
		return

	return DateArr


def Compute_Statistics(idata=[],ilats=[],ilons=[]):#,imask=[]): 
	"""
	This function computes spatial mean and RMS of a given input data.
	
	Input: 
	 	idata (lats,lons,time) input data
	 	ilats (list) latitudes in degree north
	 	ilons (list) longitudes in degree east
	Output: 
		mean, rms, min, max, area, sum  (time) 
	
	Notes:
		* mean value depends on orientation of 'lats' for whatever reason -> must be arranged from N to S
	"""
	
	import numpy as np
	from inspect import currentframe, getframeinfo

	if len(idata) == 0:
		print(f'In routine: {_currentroutine}')
		print(Compute_Statistics.__doc__)
		return

	if idata.shape[0] != len(ilats) or idata.shape[1] != len(ilons):
		frameinfo = getframeinfo(currentframe())
		print(f'ERROR! Shape of input field does not match number of lats and lons ({_currentroutine}:line {frameinfo.lineno+2}).')
		return
	
	if len(idata.shape) == 2:
		idata = idata[...,np.newaxis]
		
	nlev = idata.shape[2]	
		
	imask = idata.mask.astype(float)	
	imask = np.flip(1-imask,0)
	
	AreaPerGridPoint, TotalArea =  _CalculateAreas(ilats,ilons)
	AreaWeights = AreaPerGridPoint/TotalArea

	FIELDMEAN = np.zeros(nlev)
	FIELDRMS  = np.zeros(nlev)
	FIELDMIN  = np.zeros(nlev)
	FIELDMAX  = np.zeros(nlev)
	FIELDSUM  = np.zeros(nlev)	
	
	for i in range(nlev):
		print(f'Compute statistics for level {i+1}/{nlev}.',end='\r')
	
		## Weight meteorological field with weights
		WeightedField = np.flip(idata[:,:,i]*AreaWeights[:,:],0)

		## Compute area fraction of the masked array with respect to global area
		WeightSum = np.sum(AreaWeights[:,:]*imask[:,:,i]) # scalar

		FIELDMEAN[i] = np.sum(WeightedField[:,:]*imask[:,:,i])/WeightSum
		WeightedFieldRMS = np.flip(np.square(idata[:,:,i])*AreaWeights[:,:],0)
		FIELDRMS[i] = np.sqrt(np.sum(WeightedFieldRMS*imask[:,:,i])/WeightSum)
		FIELDMIN[i] = np.ma.min((np.float128(idata[:,:,i])))
		FIELDMAX[i] = np.ma.max((np.float128(idata[:,:,i])))
		FIELDSUM[i] = np.sum(np.flip(idata[:,:,i],0)*imask[:,:,i]*AreaPerGridPoint[:,:])

	return np.squeeze(FIELDMEAN), np.squeeze(FIELDRMS), np.squeeze(FIELDMIN), np.squeeze(FIELDMAX), TotalArea*WeightSum, np.squeeze(FIELDSUM)	


def Compute_PatternCorrelation(idata1=[],idata2=[],lats=[], detrend=False):
	"""
	This routine computes the pattern correlation between idata1 and idata2."
	
	Input:
		idata1  (2D array) input field 1
		idata2  (2D array) input field 2
		lats    (list)     latitudes of both fields
		detrend (bool)     detrending before correlation coefficients are computed
		
	Output:
		correlation coefficients (2D array)
	"""

	import numpy as np
	import scipy.signal as signal
	
	if len(idata1) == 0:
		print(f'In routine: {_currentroutine}')
		print(Compute_PatternCorrelation.__doc__)
		return
	
	if not np.ma.isMaskedArray(idata1):
		idata1 = np.ma.copy(idata1)
		idata1.mask = [False]
	
	if not np.ma.isMaskedArray(idata2):
		idata2 = np.ma.copy(idata2)
		idata2.mask = [False]

	index = np.where(idata1.mask == False) and np.where(idata2.mask == False)

	cos_phi = np.cos(lats/180*np.pi)
	idata1 = np.ma.multiply(idata1.T,cos_phi).T
	idata2 = np.ma.multiply(idata2.T,cos_phi).T

	print(f'Number of unmasked grid points: {len(index[0])}')
	
	if detrend:
		idata1_detr = signal.detrend(idata1[index])
		idata2_detr = signal.detrend(idata2[index])
		
		return np.corrcoef(idata1_detr,idata2_detr)[0,1]
	else:
		return np.corrcoef(idata1[index],idata2[index])[0,1]
		
	
	
def Compute_Tendency(idata=[],date=[]):
	"""
	This function computes the temporal tendency of a field (e.g., OHCT).
	
	Input:
		idata (lats,lons,time) 3D input array.
		date  (list)           list of datetime objects.
		
	Output:
		tendencies, with same dimensions as idata.
		
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
	
	# this part removes all indices that are masked and creates a vector of data -> for better performance
	shape = idata.shape
	mask = idata.mask[:,:,0]
	index = np.where(mask == False)
	idata = idata[index]

	odata = np.zeros_like(idata)
	
	# first timestep
	odata[:,0] = (idata[:,1] - idata[:,0])/(date[1]-date[0]).total_seconds()
	
	# last timestep
	odata[:,-1] = (idata[:,-1] - idata[:,-2])/(date[-1]-date[-2]).total_seconds()
		
	# all other timesteps
	for i in range(1,shape[-1]-1):
		odata[:,i] = (idata[:,i+1] - idata[:,i-1])/(date[i+1]-date[i-1]).total_seconds()
		
	# this part converts the vector back to an array
	odata_out = np.ma.zeros(shape)
	odata_out.mask = [True]
	odata_out.mask[index] = False
	odata_out[index] = odata
	
	return odata_out
	
	
def Compute_MeanOverTime(idata=[], DateArr=[], sdate=None, edate=None, dfmt='%Y-%m-%d', print_period=False):
	
	"""
	This function computes the temporal mean of monthly data over the period sdate-edate.
	if sdate and edate are not defined, the mean over the entire period is computed.
	
	Input: 
	 	idata   (lats,lons,time)
	 	DateArr (list of datetime objects)
	 	sdate   (str or datetime object)
	 	edate   (str or datetime object)
	 	dfmt    (str) datetime format, default is '%Y-%m-%d'
	 	
	Output:
		Temporal mean of idata with dimensions (lats,lons)
	
	Notes:
		* This function considers the length of each calender month.
		* ONLY for monthly means (not for other time resolutions, e.g., seasonal means). 
	"""
	
	import datetime
	import numpy as np


	if len(idata) == 0 or len(DateArr) == 0:
		print(f'In routine: {_currentroutine}')
		print(Compute_MeanOverTime.__doc__)
		return
		
	if sdate == None: sdate = DateArr[0]
	if edate == None: edate = DateArr[-1]

	if type(sdate) == str: sdate = datetime.datetime.strptime(sdate,dfmt)
	if type(edate) == str: edate = datetime.datetime.strptime(edate,dfmt)
	
	sindex = _FindDateIndex(sdate,DateArr)
	eindex = _FindDateIndex(edate,DateArr) + 1

	DaysPerMonth = np.array([31,28,31,30,31,30,31,31,30,31,30,31])	

	odata = np.ma.zeros(idata.shape[:2])
	cnt_days = 0
	if print_period == True: print(f'Average over {DateArr[sindex]} -- {DateArr[eindex-1]}')
	for i in range(sindex,eindex):
		if DateArr[i].year%4 == 0:
			DaysPerMonth[1] = 29
		else:
			DaysPerMonth[1] = 28 

		odata += idata[:,:,i]*DaysPerMonth[DateArr[i].month-1]
		cnt_days += DaysPerMonth[DateArr[i].month-1]

	odata = odata[:,:]/cnt_days

	return odata
	
	
def Compute_AnomalyClimatology(idata=[],DateArr=[],refdate = [None], dfmt='%Y-%m-%d',print_info=False):
	# statt 1m30s -> 9s!!
	
	"""
	This function computes the climatology and corresponding anomalies.
	
	Input:
		idata (lats,lons,time) input data  
		DateArr (list of datetime objects)
		refdate (list of str) reference period containing start and end date as strings, e.g.: ['2000-01-15','2009-12-15']
		dfmt (str) datetime object format
		
	Output:
		anomalies, climatology
		
	Notes:
		* If 'refdate' is not given, the full period is used to compute the climatology
	"""
	
	import numpy as np
	import datetime
	
	if len(idata) == 0:
		print(f'In routine: {_currentroutine}')
		print(Compute_AnomalyClimatology.__doc__)
		return

	# converts array into vector where masked indices are neglected -> better performance
	full_shape = idata.shape
	mask = idata[:,:,0].mask # maske wird nur von 1. zeitschritt übernommen!!!
	index = np.where(mask == False)

	idata = idata[index][:,np.newaxis,:]
	mask_shape = idata.shape
	
	nlat = idata.shape[0]
	nlon = idata.shape[1]
	ntime = idata.shape[2]
	anomaly = np.ma.zeros(idata.shape)
	climatology = np.ma.zeros([nlat,nlon,12])
	iorder = np.arange(12)
	
	if refdate[0] == None:
		# if no reference period is given

		if DateArr[0].month != 1:
			iorder = np.roll(iorder,-(12-DateArr[0].month+1))

		for i in range(12):
			if print_info: print(f'Compute climatology {i}/12.')
			climatology[:,:,i] = np.nanmean(idata[:,:,iorder[i]::12],axis=2)

		if print_info: print('')
		for i in range(ntime):
			if print_info: print(f'Compute anomaly {i}/{ntime}.')
			anomaly[:,:,i] = idata[:,:,i] - climatology[:,:,DateArr[i].month-1]
		
	else:
		# computes climatology with respect to a reference period
		# anomalies relative to reference climatology, but for the whole period given by 'DateArr'

		if all([type(i) == str for i in refdate]):
			ref_sdate = datetime.datetime.strptime(refdate[0],dfmt)
			ref_edate = datetime.datetime.strptime(refdate[1],dfmt)
		else:
			ref_sdate = refdata[0]
			ref_edate = refdata[1]

		si = 0
		ei = -1
		for i in range(ntime):
			if DateArr[i] == ref_sdate:
				si = i
			if DateArr[i] == ref_edate:
				ei = i+1

		ref_data = idata[:,:,si:ei]
		ref_date = DateArr[si:ei]
		print(f' -- reference date: {ref_date[0]} to {ref_date[-1]}')
				
		if ref_date[0].month != 1:
			iorder = np.roll(iorder,-(12-ref_date[0].month+1))	
		for i in range(12):
			climatology[:,:,i] = np.nanmean(ref_data[:,:,iorder[i]::12],axis=2)	
		
		for i in range(ntime):
			anomaly[:,:,i] = idata[:,:,i] - climatology[:,:,DateArr[i].month-1]

	# transformation back to array-like data
	anomaly_out = np.ma.zeros(full_shape)
	anomaly_out.mask = [True]
	anomaly_out[index] = anomaly[:,0,:]
	anomaly_out.mask[index] = False
	
	climatology_out = np.ma.zeros([full_shape[0],full_shape[1],12])
	climatology_out.mask = [True]
	climatology_out[index] = climatology[:,0,:]
	climatology_out.mask[index] = False	

	return anomaly_out, climatology_out	
		
	
def Compute_TrendSignificance(idata=[], DateArr=[], CL=0.95, trend_factor=120, ano=False, ac='ar1'):
	"""
	This function computes the trend and corresponding significance of gridded data.
	The mask of the data is considered.
	
	Input:
		idata (lats,lons,time) input data
		DateArr (list of datetime objects) time
		CL (float) confidence level, default is 0.95
		trend_factor (int) trend factor, default is 120 (decadal trend)
		ano (bool) True if input data are anomalies, default is False
		ac (str) defines how the effective sample size is computed, possible values are 'ar1' (default) and 'full'
		
	Output:
		anomalies, trend, significance (if ano==False)
		trend, significance (if ano==True)
		
	Notes:
		* Not for time-dependent masks (only the mask of the first timestep is considered)
		* if ano==False, anomalies are computed with the function Compute_AnomalyClimatology
	"""

	import numpy as np
	import scipy.stats as stats
	from sklearn import linear_model

	if len(idata) == 0 or len(DateArr) == 0:
		print(f'In routine: {_currentroutine}')
		print(Compute_TrendSignificance.__doc__)
		return

	full_shape = idata.shape
	mask = idata[:,:,0].mask # maske wird nur von 1. zeitschritt übernommen!!!
	index = np.where(mask == False)

	idata = idata[index][:,np.newaxis,:]
	shape = idata.shape

	if ano == False:
		anomalies, _ = Compute_AnomalyClimatology(idata, DateArr)
		idata_anom = anomalies.reshape([shape[0]*shape[1],shape[2]]).T
	else:
		anomalies = idata[:,:,:]
		idata_anom = idata.reshape([shape[0]*shape[1],shape[2]]).T
	
	xdata = np.arange(len(DateArr))  

	Ntrue = len(DateArr)
	alpha = 1 - CL

	model_ols =  linear_model.LinearRegression()
	model_ols.fit(xdata.reshape(-1, 1),idata_anom) 
	idata_anom_trend = model_ols.coef_[:,0] 

	SignificanceArray = np.zeros([shape[0]*shape[1]])
	for i in range(len(idata_anom_trend)):
		y_fit = model_ols.intercept_[i] + xdata[:]*idata_anom_trend[i]
		regression_residual = idata_anom[:,i] - y_fit[:]
		r = np.corrcoef(regression_residual[1:],regression_residual[:-1])[0,1]

		Neff = EffectiveSampleSize(regression_residual,ac=ac)
		residual_variance = 1/(Neff)*np.sum(regression_residual[:]**2.)
		Standard_Error = np.sqrt(residual_variance/(np.sum((xdata[:] - np.mean(xdata[:]))**2.))) 

		df = Neff - 2
		tvalue = stats.t.ppf(1. - alpha/2.,df)
		beta_0 = 0. # i.e. the null hypothesis is that x and y are uncorrelated -> no trend
		tstatistic = np.abs((idata_anom_trend[i] - beta_0)/Standard_Error)
		SignificanceArray[i] = (0 if tstatistic < tvalue else 1)
        
	trend = idata_anom_trend.reshape([shape[0],shape[1]])*trend_factor
	significance = SignificanceArray.reshape([shape[0],shape[1]])

	anomalies_out = np.ma.zeros(full_shape)
	anomalies_out.mask = [True]

	trend_out = np.ma.zeros([full_shape[0],full_shape[1]])
	trend_out.mask = [True]

	significance_out = np.ma.zeros([full_shape[0],full_shape[1]])
	significance_out.mask = [True]   

	anomalies_out[index] = anomalies[:,0,:]
	anomalies_out.mask[index] = False

	trend_out[index] = trend[:,0]
	trend_out.mask[index] = False

	significance_out[index] = significance[:,0]
	significance_out.mask[index] = False 
       
	if ano == False:
		return anomalies_out, trend_out, significance_out
	else:
		return trend_out, significance_out	
	
	
def Apply_PredefinedMask(idata=[],ilats=[],ilons=[],iregion=None,gridType=None,print_info=False):
	"""
	 Input: (idata,ilats,ilons,iregion,gridType)
	Output: idata, mask

	FOR SINGLE LEVEL DATA!
	"""

	import eccodes
	import numpy as np
	from inspect import currentframe, getframeinfo

	if len(idata) == 0 or len(ilats) == 0 or len(ilons) == 0 or iregion == None or gridType == None:
		print(f'In routine: {_currentroutine}')
		print(Apply_PredefinedMask.__doc__)
		return

	idata = np.ma.array(idata,mask=[False])
	nlat = len(ilats)
	nlon = len(ilons)
	mask1 = np.zeros([nlat,nlon])

	if gridType == 'regular_gg':
		filename = _Path_Mask+'LSM-ERA5-F'+str(nlon//4)+'.grib' #LANDSEA-MASK-F'+str(nlon//4)+'.grib')
	elif gridType == 'regular_ll':
		if nlat%2 == 0:
			filename = _Path_Mask+'LSM-ERA5-LL'+str(nlon//4)+'.grib'
		else:
			filename = _Path_Mask+'LSM-ERA5-LL'+str(nlon//4)+'-uneven.grib'
	else:
		frameinfo = getframeinfo(currentframe())
		print(f'ERROR! Land-sea mask not found ({_currentroutine}:line {frameinfo.lineno+2}).')
		return


	if print_info: print(f' -- Load mask: {filename}')
	ifile = open(filename)
	igrib = eccodes.codes_grib_new_from_file(ifile)
	mask2D = eccodes.codes_get_values(igrib).reshape([nlat,nlon])
	latOFGP_mask = eccodes.codes_get(igrib,'latitudeOfFirstGridPointInDegrees')
	eccodes.codes_release(igrib)
	ifile.close()	


	if np.float32(np.abs(round(ilats[0],3))) != np.float32(np.abs(round(latOFGP_mask,3))): 
		frameinfo = getframeinfo(currentframe())
		print(f'ERROR! Data lats ({ilats[0]:.3f}) do not match mask lats ({latOFGP_mask:.3f})({_currentroutine}:line {frameinfo.lineno+2}).')
		return		

	## per default Land maskieren
	mask1[mask2D > 0.5] = 1.

	## lakes
	mask1[_mila(nlat,30):_mila(nlat,60),  _milo(nlon,42):_milo(nlon,115)]  = 1. 	#asia
	mask1[_mila(nlat,86):_mila(nlat,105), _milo(nlon,20):_milo(nlon,39)]   = 1. 	#africa
	mask1[_mila(nlat,105):_mila(nlat,115),_milo(nlon,20):_milo(nlon,30)]   = 1. 	#africa
	mask1[_mila(nlat,117):_mila(nlat,122),_milo(nlon,135):_milo(nlon,138)] = 1. 	#australia
	mask1[_mila(nlat,23):_mila(nlat,51),  _milo(nlon,236):_milo(nlon,265)] = 1. 	#america
	mask1[_mila(nlat,40):_mila(nlat,52),  _milo(nlon,265):_milo(nlon,285)] = 1. 	#america
	
	# new 2022-10-18
	mask1[_mcla(nlat,56):_mcla(nlat,50),  _mclo(nlon,360-77.5):_mclo(nlon,360-60)] = 1. 	# single gridpoints in newfoundland (F320)
	mask1[_mcla(nlat,54):_mcla(nlat,52),  _mclo(nlon,360-60):_mclo(nlon,360-57)] = 1. 
	mask1[_mcla(nlat,52):_mcla(nlat,50),  _mclo(nlon,360-90):_mclo(nlon,360-87)] = 1.
	mask1[_mcla(nlat,14):_mcla(nlat,11),  _mclo(nlon,12):_mclo(nlon,16)] = 1. 		# africa single gridpoint F320
	mask1[_mcla(nlat,61):_mcla(nlat,58),  _mclo(nlon,11):_mclo(nlon,16)] = 1. 		# sweden single gp F320
	mask1[_mcla(nlat,8):_mcla(nlat,5),    _mclo(nlon,360-65):_mclo(nlon,360-61)] = 1. 	# south america
	mask1[_mcla(nlat,5):_mcla(nlat,2),    _mclo(nlon,360-60):_mclo(nlon,360-53)] = 1. 	# -
	mask1[_mcla(nlat,50):_mcla(nlat,45),  _mclo(nlon,360-75):_mclo(nlon,360-71)] = 1. # North america
	mask1[_mcla(nlat,69):_mcla(nlat,64),  _mclo(nlon,360-72):_mclo(nlon,360-69)] = 1. # -
	mask1[_mcla(nlat,59):_mcla(nlat,56),  _mclo(nlon,25):_mclo(nlon,29)] = 1. # estland

	if iregion in ['Europe','Asia','Australia','NorthAmerica','SouthAmerica','Antarctica','Africa']:  mask1 = 1 - mask1[:,:]
	if iregion == 'land':  mask1 = 1 - mask1[:,:]

	mask1[_DefineMask(nlat,nlon,iregion,ilats,ilons)] = 1
	idata.mask[mask1 == 1] = True
	return idata, mask1


def _mila(nlat,index):
	# mask index latitude
	# index for mask according to input INDEX 
	return int(nlat/180.*index)

def _milo(nlon,index):
	# mask index longitude
	return int(nlon/360.*index)

def _mclo(nlon,ilon):
	# index for mask according to input LONGITUDE 
	# assuming lon goes from 0 to 360
	return round(nlon/360.*ilon)


def _mcla(nlat,ilat):
	# index for mask according to input LATITUDE 
	# assuming lat goes from 90 to -90
	return round(nlat/180.*(90-ilat))	

def _macola(ilat,lats):
	import numpy as np
	# index for mask according to input LATITUDE 
	nlats = len(lats)
	idx = (np.abs(lats - ilat)).argmin()
	return idx #np.where(lats == idx) #int(nlats/180.*(lats[0]-ilat))
	

def _get_coords(input_string):
	from inspect import currentframe, getframeinfo

	if len(input_string) != 4:
		frameinfo = getframeinfo(currentframe())
		print(f'ERROR!  Wrong length of input_string ({_currentroutine}:line {frameinfo.lineno+2}).')
		exit()	
	
	lat0 = lat1 = lon0 = lon1 = None
	lat_cnt = lon_cnt = 0
	for i in range(4):
		if input_string[i][-1] in ['N','S']:
			lat_cnt += 1
			if lat_cnt > 2:
				frameinfo = getframeinfo(currentframe())
				print(f'ERROR!  Wrong coordinates in input_string ({_currentroutine}:line {frameinfo.lineno+2}).')
				exit()

			if lat0 == None:
				lat0 = float(input_string[i][:-1])
				if input_string[i][-1] == 'S': lat0 *= -1.
				continue
			if lat1 == None:
				lat1 = float(input_string[i][:-1])
				if input_string[i][-1] == 'S': lat1 *= -1.
				continue
		if input_string[i][-1] in ['E','W']:
			lon_cnt += 1
			if lon_cnt > 2:
				frameinfo = getframeinfo(currentframe())
				print(f'ERROR!  Wrong coordinates in input_string ({_currentroutine}:line {frameinfo.lineno+2}).')
				exit()
			if lon0 == None:
				lon0 = float(input_string[i][:-1])
				if input_string[i][-1] == 'W': lon0 = 360.-lon0
				continue
			if lon1 == None: 
				lon1 = float(input_string[i][:-1])
				if input_string[i][-1] == 'W': lon1 = 360.-lon1
				continue

	if lat0 == lat1 or lon0 == lon1:
		frameinfo = getframeinfo(currentframe())
		print(f'ERROR!  Wrong coordinates in input_string ({_currentroutine}:line {frameinfo.lineno+2}).')
		exit()
	if lat0 < lat1:
		tmp = lat0
		lat0 = lat1
		lat1 = tmp
	if lon1 < lon0:
		tmp = lon0
		lon0 = lon1
		lon1 = tmp	

	return [lat0,lat1,lon0,lon1]


def _FindDateIndex(sdate,DateArr, dfmt='%Y-%m-%d'):
	# sdate: datetime object
	# DateArr: Array of datetime objects

	from inspect import currentframe, getframeinfo
	from sys import exit
	
	if type(sdate) == str: sdate = datetime.datetime.strptime(sdate,dfmt)

	if not any([sdate == i for i in DateArr]):
		frameinfo = getframeinfo(currentframe())
		print(f'ERROR! Date {sdate} not found in time array ({_currentroutine}:line {frameinfo.lineno+2}).')
		exit()	

	index = 0
	for i in range(len(DateArr)):
		if sdate == DateArr[i]: index = i
	return index	
	

def _CalculateAreas(lats,lons):
	import numpy as np
	from area import area
	# returns AreaPerGridPoint, GlobalArea

	## GLOBAL AREA
	TopLeft =     [-180., 90.]
	TopRight =    [180.,  90.]
	BottomRight = [180., -90.]
	BottomLeft =  [-180, -90.]
	obj = {'type':'Polygon','coordinates':[[TopLeft,TopRight,BottomRight,BottomLeft,TopLeft]]}
	GlobalArea = area(obj)
	##

	## AREA PER GRID POINT
	StripeArea = np.zeros(len(lats))
	dlats = [abs(lats[i] - lats[i+1]) for i in range(len(lats)-1)]
	dlats = np.array(dlats)
	dlats = np.concatenate([dlats,[dlats[-1]]])
	dlon = abs(lons[1] - lons[0])
	TopLeft =  [0.,90.] 
	TopRight = [0.+dlon,90.]
	
	# Compute one meridional stripe and Copy nlon times
	for ilat in range(len(lats)):
		BottomRight = [0.+dlon, lats[ilat] - dlats[ilat]/2.]
		BottomLeft =  [0., lats[ilat] - dlats[ilat]/2.]
		if ilat == len(lats)-1:
			BottomRight = [0.+dlon, -90.]
			BottomLeft =  [0., -90.]
		obj = {'type':'Polygon','coordinates':[[TopLeft,TopRight,BottomRight,BottomLeft,TopLeft]]}
		StripeArea[ilat] = area(obj)
		TopLeft = BottomLeft[:]
		TopRight = BottomRight[:]
	
	AreaPerGridPoint = [StripeArea[:].T for i in range(len(lons))]

	return np.array(AreaPerGridPoint).T, GlobalArea
	
	
def _DefineMask(nlat,nlon,region,lats,lons):
	
	import numpy as np
	clon = 0.

	MaskArray = np.zeros([nlat,nlon])

	if region in ['Indic','indic']:
		MaskArray[:_mila(nlat,60),:]  = 1
		MaskArray[:,:_milo(nlon,22)]  = 1
		MaskArray[:,_milo(nlon,147):] = 1
		MaskArray[_mila(nlat,30):_mila(nlat,97),_milo(nlon,103):] = 1
	elif region in ['Atlantic','atlantic']:	
		MaskArray[:_mila(nlat,24),:] = 1.	# mask  arctic		
		MaskArray[_mila(nlat,82):,_milo(nlon,22):_milo(nlon,290)] = 1.	# mask south pacific
		MaskArray[_mila(nlat,58):_mila(nlat,90),_milo(nlon,30):_milo(nlon,260)]  = 1.	# masks north pacific
		MaskArray[_mila(nlat,24):_mila(nlat,58),_milo(nlon,45):_milo(nlon,270)]  = 1.
		MaskArray[_mila(nlat,73):_mila(nlat,83),_milo(nlon,260):_milo(nlon,272)] = 1.
		MaskArray[_mila(nlat,78):_mila(nlat,83),_milo(nlon,272):_milo(nlon,277)] = 1.
	elif region in ['Pacific','pacific']:	
		MaskArray[:_mila(nlat,24),:] 				= 1.	# arctic
		MaskArray[:,_milo(nlon,290):] 				= 1.	#south pacific
		MaskArray[_mila(nlat,93):,:_milo(nlon,148)] 	= 1.	# south pacific
		MaskArray[:,:_milo(nlon,99)] 				= 1.	# indic
		MaskArray[:_mila(nlat,73),_milo(nlon,260):] 	= 1.	# caribic
		MaskArray[_mila(nlat,74):_mila(nlat,82),_milo(nlon,276):] = 1.	# caribic
		MaskArray[_mila(nlat,73):_mila(nlat,75),_milo(nlon,271):_milo(nlon,290)] = 1.	# caribic 89-70 W
	elif region == 'AtlanticArctic':
		MaskArray[_mila(nlat,82):,_milo(nlon,22):_milo(nlon,290)] = 1.	# mask south pacific
		MaskArray[_mila(nlat,58):_mila(nlat,90),_milo(nlon,30):_milo(nlon,260)]  = 1.	# masks north pacific
		MaskArray[_mila(nlat,24):_mila(nlat,58),_milo(nlon,45):_milo(nlon,270)]  = 1.
		MaskArray[_mila(nlat,73):_mila(nlat,83),_milo(nlon,260):_milo(nlon,272)] = 1.
		MaskArray[_mila(nlat,78):_mila(nlat,83),_milo(nlon,272):_milo(nlon,277)] = 1.
		
	elif region == 'ARCGATE2FRAM':
		## ARCGATE DOMAIN
		# Südlich von 67N (Bering + Davis + Western GSR)
		MaskArray[_mcla(nlat,67.0):,_mclo(nlon,45):_mclo(nlon,360.-33.)] = 1.
		# Eastern GSR
		masktest = _LineSelection([60.0,-4.+360],[71.0,-40.+360],lats,lons,'S')
		MaskArray[masktest == 1] = 1
		# Alles südlich von 51N (Ärmelkanal)
		MaskArray[_mcla(nlat,51.):,:] = 1.		
		
		
		MaskArray[_mcla(nlat,79):_mcla(nlat,45),_mclo(nlon,315):_mclo(nlon,360)]  = 1.
		MaskArray[_mcla(nlat,80):_mcla(nlat,45),_mclo(nlon,0):_mclo(nlon,17)]  = 1.	
		
		# ostsee
		MaskArray[_mcla(nlat,66):_mcla(nlat,52),_mclo(nlon,14):_mclo(nlon,33)]  = 1.		
		MaskArray[_mcla(nlat,65):_mcla(nlat,55),_mclo(nlon,15.):_mclo(nlon,20.)] = 1.
		
	elif region == 'ARCGATE':
		## ARCGATE DOMAIN
		# Südlich von 67N (Bering + Davis + Western GSR)
		MaskArray[_mcla(nlat,67.0):,_mclo(nlon,45):_mclo(nlon,360.-33.)] = 1.
		# Eastern GSR
		masktest = _LineSelection([60.0,-4.+360],[71.0,-40.+360],lats,lons,'S')
		MaskArray[masktest == 1] = 1
		# Alles südlich von 51N (Ärmelkanal)
		MaskArray[_mcla(nlat,51.):,:] = 1.

	elif region == 'Atlantic2Fram':
		# Maskiere indik und pazifik
		#MaskArray[_mcla(nlat,30):,_milo(nlon,30):_milo(nlon,260)] = 1.
		
		MaskArray[_mila(nlat,82):,_milo(nlon,22):_milo(nlon,290)] = 1.	# mask south pacific
		MaskArray[_mila(nlat,58):_mila(nlat,90),_milo(nlon,30):_milo(nlon,260)]  = 1.	# masks north pacific
		MaskArray[_mila(nlat,24):_mila(nlat,58),_milo(nlon,45):_milo(nlon,270)]  = 1.
		MaskArray[_mila(nlat,73):_mila(nlat,83),_milo(nlon,260):_milo(nlon,272)] = 1.
		MaskArray[_mila(nlat,78):_mila(nlat,83),_milo(nlon,272):_milo(nlon,277)] = 1.
		
		# Mittelmeer
		MaskArray[_mcla(nlat,48):_mcla(nlat,28),_milo(nlon,0):_milo(nlon,45)] = 1. 	
		MaskArray[_mcla(nlat,40):_mcla(nlat,35),_milo(nlon,354):_milo(nlon,360)] = 1.
		# Arktis nördlich von 78.5N (Fram Strait)
		MaskArray[:_mcla(nlat,78.5),:] = 1.
		# Nordpazifik und Arktis von 30E bis 180E bis 100W
		MaskArray[:_mcla(nlat,26.5),_milo(nlon,30):_milo(nlon,260)]  = 1.	
		# Karibik Rest, westlich von Florida
		#MaskArray[_mcla(nlat,34):_mcla(nlat,25),_milo(nlon,260):_milo(nlon,278.3)] = 1.
		# Östlich von 19.5E (BSO) bis Fury and Hecla Strait
		MaskArray[_mcla(nlat,81):_mcla(nlat,66),_mclo(nlon,19.5):_mclo(nlon,180+97.)] = 1.	#_milo(nlon,30)		
		# Nördlich von ca. 67N (Devis Strait)
		masktest = _LineSelection([66.5,-74.+360],[68.5,-53.+360],lats,lons,'N')
		MaskArray[masktest == 1] = 1
		# Rest der Buffin Bay
		MaskArray[_mcla(nlat,85):_mcla(nlat,71),_mclo(nlon,360.-88.):_mclo(nlon,360.-70.)] = 1.	
	
	elif region == 'DomainN+S':
		# Maskiere alles südlich von 26.5N
		MaskArray[_mcla(nlat,26.5):,:] = 1.
		# Mittelmeer
		MaskArray[_mcla(nlat,48):_mcla(nlat,28),_milo(nlon,0):_milo(nlon,45)] = 1. 	
		MaskArray[_mcla(nlat,40):_mcla(nlat,35),_milo(nlon,354):_milo(nlon,360)] = 1.
		# Arktis nördlich von 78.5N (Fram Strait)
		MaskArray[:_mcla(nlat,78.5),:] = 1.
		# Nordpazifik und Arktis von 30E bis 180E bis 100W
		MaskArray[:_mcla(nlat,26.5),_milo(nlon,30):_milo(nlon,260)]  = 1.	
		# Karibik Rest, westlich von Florida
		MaskArray[_mcla(nlat,34):_mcla(nlat,25),_milo(nlon,260):_milo(nlon,278.3)] = 1.
		# Östlich von 19.5E (BSO) bis Fury and Hecla Strait
		MaskArray[_mcla(nlat,81):_mcla(nlat,66),_mclo(nlon,19.5):_mclo(nlon,180+97.)] = 1.	#_milo(nlon,30)		
		# Nördlich von ca. 67N (Devis Strait)
		masktest = _LineSelection([66.5,-74.+360],[68.5,-53.+360],lats,lons,'N')
		MaskArray[masktest == 1] = 1
		# Rest der Buffin Bay
		MaskArray[_mcla(nlat,85):_mcla(nlat,71),_mclo(nlon,360.-88.):_mclo(nlon,360.-70.)] = 1.
	elif region == 'DomainN':
		# Maskiere alles südlich von 26.5N
		MaskArray[_mcla(nlat,50.5):,:] = 1.
		# Arktis nördlich von 78.5N (Fram Strait)
		MaskArray[:_mcla(nlat,78.5),:] = 1.
		# Nordpazifik und Arkstis von 30E bis 180E bis 100W
		MaskArray[:_mcla(nlat,26.5),_milo(nlon,30):_milo(nlon,260)]  = 1.	
		# Östlich von 19.5E (BSO) bis Fury and Hecla Strait
		MaskArray[_mcla(nlat,81):_mcla(nlat,50),_mclo(nlon,19.5):_mclo(nlon,180+97.)] = 1.	#_milo(nlon,30)		
		# Hudson + Buffin Bay
		MaskArray[_mcla(nlat,85):_mcla(nlat,45),_mclo(nlon,360.-105.):_mclo(nlon,360.-40.)] = 1.
		# Südlich von ca. 68N (GSR)
		masktest = _LineSelection([61.0,-1.5+360],[71.0,-40.+360],lats,lons,'S')
		MaskArray[masktest == 1] = 1

		# Nordsee
		MaskArray[_mcla(nlat,61):_mcla(nlat,50),_mclo(nlon,0.):_mclo(nlon,20.)] = 1.
		MaskArray[_mcla(nlat,61):_mcla(nlat,50),_mclo(nlon,360.-1.5):_mclo(nlon,360.)] = 1.

		# Einzelne Gitterpunkte in Ostsee
		MaskArray[_mcla(nlat,65):_mcla(nlat,55),_mclo(nlon,15.):_mclo(nlon,20.)] = 1.

	elif region == 'test':
		masktest = _LineSelection([66.5,-74.+360],[68.5,-53.+360],lats,lons,'N')
		MaskArray[masktest == 1] = 1
		# Nördlich von ca. 68N (GSR)
		masktest = _LineSelection([61.0,-1.5+360],[71.0,-40.+360],lats,lons,'N')
		MaskArray[masktest == 1] = 1

	elif region in ['DomainS','domains','domainS']:
		# Maskiere alles südlich von 26.5N
		MaskArray[_mcla(nlat,26.5):,:] = 1.
		# Mittelmeer
		MaskArray[_mcla(nlat,48):_mcla(nlat,28),_milo(nlon,0):_milo(nlon,45)] = 1. 	
		MaskArray[_mcla(nlat,40):_mcla(nlat,35),_milo(nlon,354):_milo(nlon,360)] = 1.
		# Rest der Buffin Bay
		MaskArray[_mcla(nlat,85):_mcla(nlat,71),_mclo(nlon,360.-88.):_mclo(nlon,360.-70.)] = 1.
		# Karibik Rest, westlich von Florida
		MaskArray[_mcla(nlat,34):_mcla(nlat,25),_milo(nlon,260):_milo(nlon,278.3)] = 1.
		# Arktis nördlich von 78.5N (Fram Strait)
		MaskArray[:_mcla(nlat,71.5),:] = 1.
		# Nordpazifik und Arktis von 30E bis 180E bis 100W
		MaskArray[:_mcla(nlat,26.5),_milo(nlon,30):_milo(nlon,260)]  = 1.	
		# Östlich von 19.5E (BSO) bis Fury and Hecla Strait
		MaskArray[_mcla(nlat,81):_mcla(nlat,66),_mclo(nlon,19.5):_mclo(nlon,180+97.)] = 1.	#_milo(nlon,30)		
		# Nördlich von ca. 67N (Devis Strait)
		masktest = _LineSelection([66.5,-74.+360],[68.5,-53.+360],lats,lons,'N')
		MaskArray[masktest == 1] = 1
		# Nördlich von ca. 68N (GSR)
		masktest = _LineSelection([61.0,-1.5+360],[71.0,-40.+360],lats,lons,'N')
		MaskArray[masktest == 1] = 1
		# Nördlich von Nordsee
		MaskArray[_mcla(nlat,71.5):_mcla(nlat,61.),_mclo(nlon,0.):_mclo(nlon,16.)] = 1.
		MaskArray[_mcla(nlat,71.5):_mcla(nlat,61.),_mclo(nlon,358.5):_mclo(nlon,360.)] = 1.
		# Reste nördlich von Norwegen
		MaskArray[_mcla(nlat,74):_mcla(nlat,67.),_mclo(nlon,14.):_mclo(nlon,20.)] = 1.

	elif region == 'SAMBA2GSR':
		# Maskiere alles südlich von 34.5S
		MaskArray[_mcla(nlat,-34.5):,:] = 1.
		# Mittelmeer
		MaskArray[_mcla(nlat,48):_mcla(nlat,28),_milo(nlon,0):_milo(nlon,45)] = 1. 	
		MaskArray[_mcla(nlat,40):_mcla(nlat,35),_milo(nlon,354):_milo(nlon,360)] = 1.
		# Rest der Buffin Bay
		MaskArray[_mcla(nlat,85):_mcla(nlat,71),_mclo(nlon,360.-88.):_mclo(nlon,360.-70.)] = 1.
		# Karibik Rest, westlich von Florida
		#MaskArray[_mcla(nlat,34):_mcla(nlat,25),_milo(nlon,260):_milo(nlon,278.3)] = 1.
		# Arktis nördlich von 78.5N (Fram Strait)
		MaskArray[:_mcla(nlat,71.5),:] = 1.
		# Nordpazifik und Arktis von 30E bis 180E bis 100W
		MaskArray[:_mcla(nlat,26.5),_milo(nlon,30):_milo(nlon,260)]  = 1.	
		# Östlich von 19.5E (BSO) bis Fury and Hecla Strait
		MaskArray[_mcla(nlat,81):_mcla(nlat,66),_mclo(nlon,19.5):_mclo(nlon,180+1.)] = 1.	#_milo(nlon,30)		
		# Nördlich von ca. 67N (Devis Strait)
		masktest = _LineSelection([66.5,-74.+360],[68.5,-53.+360],lats,lons,'N')
		MaskArray[masktest == 1] = 1
		# Nördlich von ca. 68N (GSR)
		masktest = _LineSelection([61.0,-1.5+360],[71.0,-40.+360],lats,lons,'N')
		MaskArray[masktest == 1] = 1
		# Nördlich von Nordsee
		MaskArray[_mcla(nlat,71.5):_mcla(nlat,61.),_mclo(nlon,0.):_mclo(nlon,16.)] = 1.
		MaskArray[_mcla(nlat,71.5):_mcla(nlat,61.),_mclo(nlon,358.5):_mclo(nlon,360.)] = 1.
		# Reste nördlich von Norwegen
		MaskArray[_mcla(nlat,74):_mcla(nlat,67.),_mclo(nlon,14.):_mclo(nlon,20.)] = 1.


		# from Atlantic region -> not checked for redundant masking
		MaskArray[:_mila(nlat,24),:] = 1.	# mask  arctic		
		MaskArray[_mila(nlat,82):,_milo(nlon,22):_milo(nlon,290)] = 1.	# mask south pacific
		MaskArray[_mila(nlat,58):_mila(nlat,90),_milo(nlon,30):_milo(nlon,260)]  = 1.	# masks north pacific
#		MaskArray[_mila(nlat,24):_mila(nlat,58),_milo(nlon,45):_milo(nlon,270)]  = 1.
		MaskArray[_mila(nlat,73):_mila(nlat,83),_milo(nlon,260):_milo(nlon,272)] = 1.
		MaskArray[_mila(nlat,78):_mila(nlat,83),_milo(nlon,272):_milo(nlon,277)] = 1.

		# single dots in SA
		MaskArray[_mcla(nlat,-19):_mcla(nlat,-24),_mclo(nlon,289):_mclo(nlon,294)] = 1.
		MaskArray[_mcla(nlat,13):_mcla(nlat,10),_mclo(nlon,270):_mclo(nlon,275)] = 1.
		MaskArray[_mcla(nlat,8.5):_mcla(nlat,7),_mclo(nlon,278):_mclo(nlon,283)] = 1.
		MaskArray[_mcla(nlat,9.5):_mcla(nlat,7.5),_mclo(nlon,286.5):_mclo(nlon,289.5)] = 1.
		MaskArray[_mcla(nlat,-29):_mcla(nlat,-32),_mclo(nlon,296):_mclo(nlon,300)] = 1.

	elif region == 'SAMBA2Bering':
		# Maskiere alles südlich von 34.5S
		MaskArray[_mcla(nlat,-34.5):,:] = 1.
		# Mittelmeer
		MaskArray[_mcla(nlat,48):_mcla(nlat,28),_milo(nlon,0):_milo(nlon,45)] = 1. 	
		MaskArray[_mcla(nlat,40):_mcla(nlat,35),_milo(nlon,354):_milo(nlon,360)] = 1.

		# from Atlantic region -> not checked for redundant masking
#		MaskArray[:_mila(nlat,24),:] = 1.	# mask  arctic		
		MaskArray[_mcla(nlat,66.):_mcla(nlat,30.),_mclo(nlon,110):_mclo(nlon,250)] = 1.
		MaskArray[_mila(nlat,82):,_milo(nlon,22):_milo(nlon,290)] = 1.	# mask south pacific
		MaskArray[_mila(nlat,58):_mila(nlat,90),_milo(nlon,30):_milo(nlon,260)]  = 1.	# masks north pacific
#		MaskArray[_mila(nlat,24):_mila(nlat,58),_milo(nlon,45):_milo(nlon,270)]  = 1.
		MaskArray[_mila(nlat,73):_mila(nlat,83),_milo(nlon,260):_milo(nlon,272)] = 1.
		MaskArray[_mila(nlat,78):_mila(nlat,83),_milo(nlon,272):_milo(nlon,277)] = 1.

		# single dots in SA
		MaskArray[_mcla(nlat,-19):_mcla(nlat,-24),_mclo(nlon,289):_mclo(nlon,294)] = 1.
		MaskArray[_mcla(nlat,13):_mcla(nlat,10),_mclo(nlon,270):_mclo(nlon,275)] = 1.
		MaskArray[_mcla(nlat,8.5):_mcla(nlat,7),_mclo(nlon,278):_mclo(nlon,283)] = 1.
		MaskArray[_mcla(nlat,9.5):_mcla(nlat,7.5),_mclo(nlon,286.5):_mclo(nlon,289.5)] = 1.
		MaskArray[_mcla(nlat,-29):_mcla(nlat,-32),_mclo(nlon,296):_mclo(nlon,300)] = 1.


	elif region in ['gulf','Gulf']:
		# Maskiere alles südlich von 26.5N
		MaskArray[_mcla(nlat,26.5):,:] = 1.
		# Mittelmeer
		MaskArray[_mcla(nlat,48):_mcla(nlat,28),_milo(nlon,0):_milo(nlon,45)] = 1. 	
		MaskArray[_mcla(nlat,40):_mcla(nlat,35),_milo(nlon,354):_milo(nlon,360)] = 1.
		# Rest der Buffin Bay
		MaskArray[_mcla(nlat,85):_mcla(nlat,71),_mclo(nlon,360.-88.):_mclo(nlon,360.-70.)] = 1.
		# Karibik Rest, westlich von Florida
		MaskArray[_mcla(nlat,34):_mcla(nlat,25),_milo(nlon,260):_milo(nlon,278.3)] = 1.
		# Arktis nördlich von 78.5N (Fram Strait)
		MaskArray[:_mcla(nlat,71.5),:] = 1.
		# Nordpazifik und Arktis von 30E bis 180E bis 100W
		MaskArray[:_mcla(nlat,26.5),_milo(nlon,30):_milo(nlon,260)]  = 1.	
		# Östlich von 19.5E (BSO) bis Fury and Hecla Strait
		MaskArray[_mcla(nlat,81):_mcla(nlat,66),_mclo(nlon,19.5):_mclo(nlon,180+97.)] = 1.	#_milo(nlon,30)		
		# Devis Strait
		MaskArray[:_mcla(nlat,55.),:] = 1.
		#
		masktest = _LineSelection([24,-56+360],[56.,-15.+360],lats,lons,'E')
		MaskArray[masktest == 1] = 1
		MaskArray[_mcla(nlat,60):_mcla(nlat,45),_mclo(nlon,0):_mclo(nlon,25)] = 1.
		
		MaskArray[_mcla(nlat,58):_mcla(nlat,48),_mclo(nlon,-85+360):_mclo(nlon,-75+360)] = 1.

	elif region == 'Africa':
		# Alles nördlich von 38N
		MaskArray[:_mcla(nlat,38.),:] = 1. 
		# alles südlich von 60S
		MaskArray[_mcla(nlat,-60.):,:] = 1. 
		# Alles östlich von 75E bis 30W
		MaskArray[:,_milo(nlon,75):_milo(nlon,360.-30.)] = 1. 
		# Teile Europas
		MaskArray[_mcla(nlat,38.):_mcla(nlat,36.),_mclo(nlon,12):] = 1. 
		# Saudi Arabien
		masktest = _LineSelection([12.,44.],[34.,30.],lats,lons,'N')
		MaskArray[masktest == 1] = 1
		MaskArray[_mcla(nlat,37.):_mcla(nlat,12.),_mclo(nlon,44.):_mclo(nlon,75.)] = 1. 
		
	elif region == 'Antarctica':
		# Alles nördlich von 60S
		MaskArray[:_mcla(nlat,-60.),:] = 1. 
				
	elif region == 'SouthAmerica' or region == 'SA':
		# Alles östlich und westlich von SA
		MaskArray[:,:_mclo(nlon,360.-83.)] = 1.
		MaskArray[:,_mclo(nlon,360.-30.):] = 1.		
		# Alles nördlich von 14N
		MaskArray[:_mcla(nlat,16.),:] = 1.
		# Alles südlich von 60S
		MaskArray[_mcla(nlat,-60.):,:] = 1.
		# Panama
		MaskArray[_mcla(nlat,10.):_mcla(nlat,5.),_mclo(nlon,360.-85.):_mclo(nlon,360.-78.)] = 1.
		
	elif region == 'NorthAmerica' or region == 'NA':
		MaskArray[_mcla(nlat,15.):_mcla(nlat,5.),_mclo(nlon,360.-78.):] = 1.
		# Alles südlich von 5N
		MaskArray[_mcla(nlat,5.):,:] = 1.
		# Island Afrika bis 0W
		MaskArray[_mcla(nlat,68.):_mcla(nlat,15.),_mclo(nlon,360.-25.):] = 1.
		#
		MaskArray[:_mcla(nlat,5.),:_mclo(nlon,360.-169.)] = 1.
		
	elif region == 'Australia':
		# südlich von 60S
		MaskArray[_mcla(nlat,-60.):,:] = 1.
		# nördlich von AQ
		MaskArray[:_mcla(nlat,0.),:] = 1.
		# Africa
		MaskArray[_mcla(nlat,0.):_mcla(nlat,-60.),:_mclo(nlon,75.)] = 1.
		# SA
		MaskArray[_mcla(nlat,0.):_mcla(nlat,-60.),_mclo(nlon,360.-90.):] = 1.
		# Philippines
		MaskArray[_mcla(nlat,0.):_mcla(nlat,-10.),_mclo(nlon,90.):_mclo(nlon,132.)] = 1.
		
	elif region == 'Asia':
		# östlich von 169W
		MaskArray[:,_mclo(nlon,360.-169.):] = 1.
		# südlich von 10S
		MaskArray[_mcla(nlat,-10.):,:] = 1.
		# Papa Neuguinea
		MaskArray[_mcla(nlat,5.):_mcla(nlat,-10.),_mclo(nlon,130.):_mclo(nlon,180.)] = 1.
		# Africa
		masktest = _LineSelection([10.,44.],[34.,29.],lats,lons,'W')
		MaskArray[masktest == 1] = 1
		MaskArray[_mcla(nlat,12.):_mcla(nlat,-10.),:_mclo(nlon,55.)] = 1.
		# Europe
		MaskArray[:_mcla(nlat,34.),:_mclo(nlon,26.)] = 1.
		masktest = _LineSelection([40.,23.],[45.,35.],lats,lons,'N')
		MaskArray[masktest == 1] = 1
		masktest = _LineSelection([40.,51.],[45.,35.],lats,lons,'N')
		MaskArray[masktest == 1] = 1
		masktest = _LineSelection([40.,51.],[76.,69.],lats,lons,'N')
		MaskArray[masktest == 1] = 1
		
	elif region == 'Europe':
		# Asia
		masktest = _LineSelection([40.,23.],[45.,35.],lats,lons,'S')
		MaskArray[masktest == 1] = 1
		masktest = _LineSelection([40.,51.],[45.,35.],lats,lons,'S')
		MaskArray[masktest == 1] = 1
		masktest = _LineSelection([40.,51.],[76.,69.],lats,lons,'S')
		MaskArray[masktest == 1] = 1
		# östlich von 69E
		MaskArray[:,_mclo(nlon,69.):_mclo(nlon,360.-30.)] = 1.
		# südlich von 36N
		MaskArray[_mcla(nlat,36.):,:] = 1.
		# greenland
		MaskArray[:_mcla(nlat,67.):,_mclo(nlon,360.-30.):_mclo(nlon,360.-10.)] = 1.
		# Rest von Afrika (Tunesien, Algerien)
		MaskArray[_mcla(nlat,38.):_mcla(nlat,35.),:_mclo(nlon,12.)] = 1.
		
	elif region == 'SAMBA2RAPID':
		# maskiere alles nördlich von 26N
		MaskArray[:_mcla(nlat,26.5),_mclo(nlon,360-81):] = 1.
		MaskArray[:_mcla(nlat,35),:] = 1.

		# Maskiere alles südlich von 34.5S
		MaskArray[_mcla(nlat,-34.5):,:] = 1.
		# Mittelmeer
		MaskArray[_mcla(nlat,48):_mcla(nlat,28),_milo(nlon,0):_milo(nlon,45)] = 1. 	
		MaskArray[_mcla(nlat,40):_mcla(nlat,35),_milo(nlon,354):_milo(nlon,360)] = 1.

		# from Atlantic region -> not checked for redundant masking
#		MaskArray[:_mila(nlat,24),:] = 1.	# mask  arctic		
		MaskArray[_mcla(nlat,66.):_mcla(nlat,30.),_mclo(nlon,110):_mclo(nlon,250)] = 1.
		MaskArray[_mila(nlat,82):,_milo(nlon,22):_milo(nlon,290)] = 1.	# mask south pacific
		MaskArray[_mila(nlat,58):_mila(nlat,90),_milo(nlon,30):_milo(nlon,260)]  = 1.	# masks north pacific
#		MaskArray[_mila(nlat,24):_mila(nlat,58),_milo(nlon,45):_milo(nlon,270)]  = 1.
		MaskArray[_mila(nlat,73):_mila(nlat,83),_milo(nlon,260):_milo(nlon,272)] = 1.
		MaskArray[_mila(nlat,78):_mila(nlat,83),_milo(nlon,272):_milo(nlon,277)] = 1.

		# single dots in SA
		MaskArray[_mcla(nlat,-19):_mcla(nlat,-24),_mclo(nlon,289):_mclo(nlon,294)] = 1.
		MaskArray[_mcla(nlat,13):_mcla(nlat,10),_mclo(nlon,270):_mclo(nlon,275)] = 1.
		MaskArray[_mcla(nlat,8.5):_mcla(nlat,7),_mclo(nlon,278):_mclo(nlon,283)] = 1.
		MaskArray[_mcla(nlat,9.5):_mcla(nlat,7.5),_mclo(nlon,286.5):_mclo(nlon,289.5)] = 1.
		MaskArray[_mcla(nlat,-29):_mcla(nlat,-32),_mclo(nlon,296):_mclo(nlon,300)] = 1.	

	return MaskArray == 1	
	
	
def _LineSelection(sarr,earr,lats,lons,MaskRegion):
	# sarr: start point [lat,lon]
	# earr: end point [lat,lon]
	# lats: array with data latitudes
	# lons: array with data longitudes
	# MaskRegion: N/E/S/W/R/L

	import numpy as np
	slat = sarr[0]
	slon = sarr[1]
	elat = earr[0]
	elon = earr[1]

	if slat > elat:
		print('ERROR! Startpunkt nördlich von Endpunkt!')
		exit()

	if MaskRegion not in ['N','S','E','W']:
		print('ERROR! Falsche Mask Region!')
		exit()

	lat_line = []
	lon_line = []
	idx_lat_NN = []
	idx_lon_NN = []
	sign_lat = int((elat-slat)/np.abs(elat-slat))
	sign_lon = int((elon-slon)/np.abs(elon-slon))

	## SELEKTIERE START GITTERPUNKT
	silat = 0	# latitude array index first gridpoint
	silon = 0	# longitude array index first gridpoint
	# wähle latitude des Start-GP aus 
	for i in range(len(lats)):
		if lats[i] < slat: 
			silat = i
			break
	sglat = lats[silat]
	if np.abs(lats[silat] - slat) > np.abs(lats[silat-1] - slat):
		sglat = lats[silat-1]

	# wähle longtiude des Start-GP aus 
	for i in range(len(lons)):
		if lons[i] > slon: 
			silon = i
			break
	sglon = lons[silon]	
	if np.abs(lons[silon] - slon) > np.abs(lons[silon-1] - slon): 
		sglon = lons[silon-1]

	## SELEKTIERE END GITTERPUNKT
	eilat = 0	# latitude array index first gridpoint
	eilon = 0	# longitude array index first gridpoint
	# wähle latitude des Start-GP aus 
	for i in range(len(lats)):
		if lats[i] < elat: 
			eilat = i
			break
	eglat = lats[eilat]
	if np.abs(lats[eilat] - elat) > np.abs(lats[eilat-1] - elat): eglat = lats[eilat-1]

	# wähle longtiude des Start-GP aus 
	for i in range(len(lons)):
		if lons[i] > elon: 
			eilon = i
			break
	eglon = lons[eilon]
	if np.abs(lons[eilon] - elon) > np.abs(lons[eilon-1] - elon): eglon = lons[eilon-1]

	nlat = np.abs(eilat-silat)+1
	nlon = np.abs(eilon-silon)+1

	if nlat >= nlon:
		slope = (slon-elon)/(slat-elat)
		lat_line += [lats[silat-1]]
		lon_line += [slon + np.abs(lats[silat-1]-slat)*slope]

		for i in range(nlat-2):
			lat_line += [lats[silat-(i+2)]]
			lon_line += [lon_line[-1] + np.abs(lat_line[-1]-lat_line[-2])*slope]
			idx_lat_NN += [int(silat-(i+2))]
		idx_lat_NN += [int(silat-(nlat))]
		lon_NN = [sglon]
		lat_NN = [sglat]

		for i in range(len(lon_line)):
			loninn = (np.abs(lons - lon_line[i])).argmin()
			lon_NN += [lons[loninn]]
			lat_NN += [lat_line[i]]
			idx_lon_NN += [int(loninn)]			
			# siehe https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array

	elif nlat < nlon:
		slope = np.abs((slat-elat)/(slon-elon))

		c1 = 0
		if sign_lon < 0.: c1 = -1

		lon_line += [lons[silon+c1]]	
		lat_line += [slat + np.abs(lons[silon]-slon)*slope]
			
		for i in range(nlon-1):
			lon_line += [lons[silon+(i+1)*sign_lon]]
			lat_line += [lat_line[-1] + np.abs(lon_line[-1]-lon_line[-2])*slope]
			idx_lon_NN += [int(silon+(i+1)*sign_lon)]

		idx_lon_NN += [int(silon+(nlon)*sign_lon)]
		lon_NN = [sglon]
		lat_NN = [sglat]
		for i in range(len(lat_line)):
			latinn = (np.abs(lats - lat_line[i])).argmin()
			lat_NN += [lats[latinn]]
			lon_NN += [lon_line[i]]
			idx_lat_NN += [int(latinn)]	
	else:
		print('ERROR!')
		exit()

	check_figure = 0
	if check_figure == 1:
		plt.figure()
		plt.plot(sglon,sglat,'ro')
		for i in range(len(lon_line)):
			plt.plot(lon_line[i],lat_line[i],'ko',alpha=0.2)
		for i in range(len(lon_NN)):
			plt.plot(lon_NN[i],lat_NN[i],'k+')
		plt.plot([slon,elon],[slat,elat],'grey',alpha=0.2)
		plt.plot(slon,slat,'rx',alpha=0.2,label='start')
		plt.plot(elon,elat,'gx',alpha=0.2,label='end')
		plt.xlim(np.min([slon,elon])-1,np.max([slon,elon])+1)
		plt.ylim(np.min([slat,elat])-1,np.max([slat,elat])+1)
		plt.legend(fontsize=8,bbox_to_anchor=(0.5,0.95))
		plt.show()
		exit()

	MaskArray = np.zeros([len(lats),len(lons)])
	for i in range(len(idx_lat_NN)):
			if MaskRegion == 'N':
				MaskArray[:idx_lat_NN[i]+1,idx_lon_NN[i]] = 1
			elif MaskRegion == 'S':	
				MaskArray[idx_lat_NN[i]+1:,idx_lon_NN[i]] = 1
			elif MaskRegion == 'E':	
				MaskArray[idx_lat_NN[i],idx_lon_NN[i]:] = 1
			elif MaskRegion == 'W':	
				MaskArray[idx_lat_NN[i],:idx_lon_NN[i]+1] = 1

	return MaskArray	
	
