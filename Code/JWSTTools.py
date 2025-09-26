#Some convenience functions for working with JWST data
#Written by R. Levy (rlevy.astro@gmail.com)
#updated for SASP2025 on 2025-07-22

class jwst_tools:
	import astropy.units as u

	def NIRCam_Filter_Info(this_filter=None):
		#info from Tables 2 and 3 of https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters#NIRCamFilters-Filterlists
		#and Table 1 of https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-point-spread-functions#NIRCamPointSpreadFunctions-NIRCamPSFs
		f140m = {'Filter': 'F140M',
				 'Pivot_microns': 1.404,
				 'BW_microns': 0.142,
				 'Effective_Response': 0.434,
				 'Half_Power_Wavelengths': [1.331, 1.479],
				 'Alt_Filter':'clear',
				 'Pixel_Scale_arcsec': 0.021,
				 'PSF_FWHM_arcsec':0.048}
		f164n = {'Filter': 'F164N',
				 'Pivot_microns': 1.644,
				 'BW_microns': 0.020,
				 'Effective_Response': 0.385,
				 'Half_Power_Wavelengths': [1.635, 1.653],
				 'Alt_Filter':'f150w2',
				 'Pixel_Scale_arcsec': 0.021,
				 'PSF_FWHM_arcsec':0.056}
		f212n = {'Filter': 'F212N',
				 'Pivot_microns': 2.120,
				 'BW_microns': 0.027,
				 'Effective_Response': 0.420,
				 'Half_Power_Wavelengths': [2.109, 2.134],
				 'Alt_Filter':'clear',
				 'Pixel_Scale_arcsec': 0.021,
				 'PSF_FWHM_arcsec':0.072}
		f250m = {'Filter': 'F250M',
				 'Pivot_microns': 2.503,
				 'BW_microns': 0.181,
				 'Effective_Response': 0.370,
				 'Half_Power_Wavelengths': [2.412, 2.595],
				 'Alt_Filter':'clear',
				 'Pixel_Scale_arcsec': 0.042,
				 'PSF_FWHM_arcsec':0.085}
		f335m = {'Filter': 'F335M',
				 'Pivot_microns': 3.365,
				 'BW_microns': 0.347,
				 'Effective_Response': 0.480,
				 'Half_Power_Wavelengths': [3.177, 3.537],
				 'Alt_Filter':'clear',
				 'Pixel_Scale_arcsec': 0.042,
				 'PSF_FWHM_arcsec':0.111}
		f360m = {'Filter': 'F360M',
				 'Pivot_microns': 3.621,
				 'BW_microns': 0.372,
				 'Effective_Response': 0.515,
				 'Half_Power_Wavelengths': [3.426, 3.814],
				 'Alt_Filter':'clear',
				 'Pixel_Scale_arcsec': 0.042,
				 'PSF_FWHM_arcsec':0.120}
		filters = {'F140M':f140m,'F164N':f164n,'F212N':f212n,'F250M':f250m,'F335M':f335m,'F360M':f360m}

		if this_filter == None:
			f = {}
		else:
			try:
				f = filters[this_filter.upper()]
			except KeyError:
				print('Invalid filter name.')
				f = {}

		return f,filters

	def MIRI_Filter_Info(this_filter=None):
		import inspect
		# get the frame object of the function
		frame = inspect.currentframe()
		func_name = frame.f_code.co_name	
		#info from Table 1 on https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-instrumentation/miri-filters-and-dispersers#MIRIFiltersandDispersers-Imagingfilters
		#and https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-observing-modes/miri-imaging
		f560w = {'Filter': 'F560W',
				 'Pivot_microns': 5.635,
				 'BW_microns': 1.0,
				 'Effective_Response': 0.245,
				 'Half_Power_Wavelengths': [5.054, 6.171],
				 'Alt_Filter':'clear',
				 'Pixel_Scale_arcsec': 0.11,
				 'PSF_FWHM_arcsec':0.207}
		f770w = {'Filter': 'F770W',
				 'Pivot_microns': 7.639,
				 'BW_microns': 1.95,
				 'Effective_Response': 0.355,
				 'Half_Power_Wavelengths': [6.581, 8.687],
				 'Alt_Filter':'clear',
				 'Pixel_Scale_arcsec': 0.11,
				 'PSF_FWHM_arcsec':0.269}
		f1130w = {'Filter': 'F1130W',
				 'Pivot_microns': 11.309,
				 'BW_microns': 0.730,
				 'Effective_Response': 0.412,
				 'Half_Power_Wavelengths': [10.953, 11.667],
				 'Alt_Filter':'clear',
				 'Pixel_Scale_arcsec': 0.11,
				 'PSF_FWHM_arcsec':0.375}
		filters = {'F560W':f560w,'F770W':f770w,'F1130W':f1130w}	

		if this_filter == None:
			f = {}
		else:
			try:
				f = filters[this_filter.upper()]
			except KeyError:
				print('Invalid filter name in %s.%s' %(__name__,func_name))
				f = {}
				exit(1)

		return f,filters

	def getNIRCamFilterNames(this_filter=None):
		from JWSTTools import jwst_tools
		one_filter,all_filters = jwst_tools.NIRCam_Filter_Info(this_filter=this_filter)
		if this_filter==None:
			filter_names = ['%s' %el.lower() for el in all_filters]
			filter_names_full = ['%s-%s' %(all_filters[el]['Alt_Filter'],el.lower()) for el in all_filters]
		else:
			filter_names = one_filter['Filter'].lower()
			filter_names_full = one_filter['Alt_Filter']+'-'+one_filter['Filter'].lower()
		return filter_names,filter_names_full

	def getMIRIFilterNames(this_filter=None):
		from JWSTTools import jwst_tools
		one_filter,all_filters = jwst_tools.MIRI_Filter_Info(this_filter=this_filter)
		if this_filter==None:
			filter_names = ['%s' %el.lower() for el in all_filters]
			filter_names_full = ['%s-%s' %(all_filters[el]['Alt_Filter'],el.lower()) for el in all_filters]
		else:
			filter_names = one_filter['Filter'].lower()
			filter_names_full = one_filter['Alt_Filter']+'-'+one_filter['Filter'].lower()
		return filter_names,filter_names_full

	def getFileNames(this_instrument=None,this_filter=None):
		from os import listdir
		from os.path import isfile, join
		from JWSTTools import jwst_tools
		import inspect
		# get the frame object of the function
		frame = inspect.currentframe()
		func_name = frame.f_code.co_name	

		data_dir = '../Data/JWST/'
		file_names = [f for f in listdir(data_dir) if isfile(join(data_dir, f))]
		#remove any non-fits files (e.g., .DS_store)
		file_names = [x for x in file_names if x.endswith('.fits')]
		#now order by increasing filter wavelength
		file_names.sort(key=lambda x: x.split('_')[3])
		#if instrument specified, only return those file names
		if this_instrument!=None:
			file_names = list(filter(lambda x: this_instrument.lower() in x, file_names))
			if not file_names:
				print('Invalid instrument name in %s.%s' %(__name__,func_name))
		#if filter specified, return only that file name
		if this_filter!=None:
			file_names = list(filter(lambda x: this_filter.lower() in x, file_names))
			if not file_names:
				print('Invalid filter name in %s.%s' %(__name__,func_name))
				exit(1)
		return file_names

	def mask_outside_footprint(ref_wcs,og_wcs,ref_shape,og_shape):
		#THIS ASSUMES THE PIXEL SCALE DOESN'T CHANGE
		import numpy as np
		import astropy.units as u
		from regions import PixCoord,RectanglePixelRegion
		#find rotation angle between orignal and new wcs
		PC = og_wcs.wcs.pc
		PA = np.degrees(np.arctan2(PC[0,1],PC[1,1]))
		#make a mask centered on the new wcs center with the size and rotation for the og wcs
		cen_pix = PixCoord(ref_wcs.wcs.crpix[0],ref_wcs.wcs.crpix[1])
		reg = RectanglePixelRegion(center=cen_pix,width=og_shape[1],height=og_shape[0],angle=PA*u.deg)
		mask = reg.to_mask()
		mask_im = mask.multiply(np.ones(ref_shape))
		mask_im[mask_im == 1.]=np.nan
		return mask_im

	def MJysr_to_MJypix(I_MJysr,hdr):
		import numpy as np
		x_pix = hdr['CDELT1'] #deg
		y_pix = hdr['CDELT2'] #deg

		pix2sr = x_pix*y_pix*np.pi/180

		return I_MJysr*pix2sr

	def Jybeam_to_Jypix(I_Jybeam,hdr):
		import numpy as np
		beam = np.pi*hdr['BMAJ']*hdr['BMIN']/(4*np.log(2))

		x_pix = hdr['CDELT1'] #deg
		y_pix = hdr['CDELT2'] #deg
		Apix = x_pix*y_pix

		return I_Jybeam*beam*Apix

	def app2abs_mag(m,dist):
		import numpy as np
		import astropy.units as u
		M = m-5*np.log10(dist.to(u.pc).value)+5
		return M
	
	def abs2app_mag(M,dist):
		import numpy as np
		import astropy.units as u
		m = M-5+5*np.log10(dist.to(u.pc).value)
		return m

	def flux2abmag(filter_name,flux):#,app_area_sr,data_version):
		import numpy as np
		from astropy.io import fits
		from JWSTTools import jwst_tools
		#convert flux in mJy to AB mag
		#https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-absolute-flux-calibration-and-zeropoints#NIRCamAbsoluteFluxCalibrationandZeropoints-ABmagnitudes
		# sb = f/npix*1E9 #MJy/pix
		# mag_ab = -6.10 - 2.5*np.log10(sb)
		# file_names = jwst_tools.getFileNames(data_version)
		# this_im = [f for f in file_names if filter_name in f][0]
		# fits_hdr = fits.open('../Data/NIRCam_Center/'+data_version+'/'+this_im)[1].header
		# pixar_sr = fits_hdr['PIXAR_SR']

		filter_name = filter_name.lower()
		nircam_filters,_ = jwst_tools.getNIRCamFilterNames()
		miri_filters,_ = jwst_tools.getMIRIFilterNames()
		if filter_name in nircam_filters:
			filter_info,_ = jwst_tools.NIRCam_Filter_Info(this_filter=filter_name)
		elif filter_name in miri_filters:
			filter_info,_ = jwst_tools.MIRI_Filter_Info(this_filter=filter_name)
		else:
			print('Invalid filter name.')
			filter_info = None

		if filter_info != None:
			pix = filter_info['Pixel_Scale_arcsec']/206265
			pixar_sr = pix**2
			zp_ab = -6.10-2.5*np.log10(pixar_sr)
			app_area_sr = pixar_sr

			sb = flux/1E9/app_area_sr #MJy/sr
			mag_ab = zp_ab - 2.5*np.log10(sb)
		else:
			mag_ab = None

		return mag_ab

	def abmag2flux(filter_name,mag_ab):#,app_area_sr,data_version):
		import numpy as np
		from astropy.io import fits
		from JWSTTools import jwst_tools
		#convert flux in mJy to AB mag
		#https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-absolute-flux-calibration-and-zeropoints#NIRCamAbsoluteFluxCalibrationandZeropoints-ABmagnitudes
		# sb = f/npix*1E9 #MJy/pix
		# mag_ab = -6.10 - 2.5*np.log10(sb)
		# file_names = jwst_tools.getFileNames(data_version)
		# this_im = [f for f in file_names if filter_name in f][0]
		# fits_hdr = fits.open('../Data/NIRCam_Center/'+data_version+'/'+this_im)[1].header
		# pixar_sr = fits_hdr['PIXAR_SR']

		filter_name = filter_name.lower()
		nircam_filters,_ = jwst_tools.getNIRCamFilterNames()
		miri_filters,_ = jwst_tools.getMIRIFilterNames()
		if filter_name in nircam_filters:
			filter_info,_ = jwst_tools.NIRCam_Filter_Info(this_filter=filter_name)
		elif filter_name in miri_filters:
			filter_info,_ = jwst_tools.MIRI_Filter_Info(this_filter=filter_name)
		else:
			print('Invalid filter name.')
			filter_info = None

		if filter_info != None:
			pix = filter_info['Pixel_Scale_arcsec']/206265
			pixar_sr = pix**2
			zp_ab = -6.10-2.5*np.log10(pixar_sr)
			app_area_sr = pixar_sr

			sb = 10**(-(mag_ab-zp_ab)/2.5) #MJy/sr
			flux = sb*1E9*app_area_sr #mJy	
		else:	
			sb = None
			flux = None

		return sb,flux
	
	def Wnm_to_mJy(Wnm,wl,dist):
		import numpy as np
		#convert W/nm to mJy
		#wl in microns
		#dist in Mpc
		d_m = dist*3.086E22 #m
		c = 2.9979E17 #nm/s
		wl_nm = wl/1E3
		return 1E29/(4*np.pi)*Wnm*wl_nm**2*d_m**-2

	def intersectionArea(x1,y1,r1,x2,y2,r2,wcs):
		import numpy as np
		if (type(x1)!=np.float64) and (x1.unit == 'deg'):
			#convert everything to pixel units
			x1,y1 = wcs.all_world2pix(x1,y1,1,ra_dec_order=True)
			x2,y2 = wcs.all_world2pix(x2,y2,1,ra_dec_order=True)
			arcsec2pix = wcs.wcs.cdelt[1]*3600
			r1 = r1.value/arcsec2pix
			r2 = r2.value/arcsec2pix

		#https://www.geeksforgeeks.org/area-of-intersection-of-two-circles/
		d = np.sqrt(((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1)))
		frac = np.zeros_like(d)

		for i in range(len(d)):
			if (d[i] > r1+r2[i]) and (d[i]>0.):
				ans = 0.
				a = 1
			elif (d[i] <= (r1-r2[i])) and (r1 >= r2[i]):
				ans = np.pi*r2[i]**2
				a = np.pi*r2[i]**2
			elif (d[i] <= (r2[i]-r1)) and (r2[i] >= r1):
				ans = np.pi*r1**2
				a = np.pi*r1**2
			else:
				alpha = np.arccos(((r1**2)+(d[i]**2)-(r2[i]**2))/(2*r1*d[i]))*2
				beta = np.arccos(((r2[i]**2)+(d[i]**2)-(r1**2))/(2*r2[i]*d[i]))*2
				a1 = (0.5*beta*r2[i]**2)-(0.5*r2[i]**2*np.sin(beta))
				a2 = (0.5*alpha*r1**2)-(0.5*r1**2*np.sin(alpha))
				ans = a1+a2
				#now return as a fraction of the smaller region area
				if r2[i] > r1:
					a = np.pi*r1**2
				elif r2[i] < r1:
					a = np.pi*r2[i]**2
				else:
					a = np.pi*r1**2

			frac[i] = ans/a
			# #only flag the smaller circle for removal!
			# if r2[i] <= r1:
			# 	frac[i] = ans/a
			# else:
			# 	frac[i] = -1

		return frac

	def cat_from_reg(reg,ra_shift=0*u.deg,dec_shift=0*u.deg):
		#convert from a region  file to a pandas dataframe catalog
		import pandas as pd
		from regions import Regions
		import numpy as np
		import astropy.units as u

		cat = pd.DataFrame(columns=['ID','RA','Dec','Radius_arcsec'])
		if ('label' in reg[0].meta):
			cat['ID'] = [reg[i].meta['label'] for i in range(len(reg))]
		elif ('text' in reg[0].meta):
			cat['ID'] = [reg[i].meta['text'] for i in range(len(reg))]
		else:
			cat['ID'] = np.arange(0,len(reg),1)
		cat['RA'] = [reg[i].center.ra.value+ra_shift.value for i in range(len(reg))]
		cat['Dec'] = [reg[i].center.dec.value+dec_shift.value for i in range(len(reg))]
		cat['Radius_arcsec'] = [reg[i].radius.to(u.arcsec).value for i in range(len(reg))]
		return cat

	def reg_from_cat(cat,color="green",fontsize=0,label=False):
		import numpy as np
		import pandas as pd
		from astropy.coordinates import SkyCoord
		import astropy.units as u
		from regions import CircleSkyRegion,Regions
		
		regions = []
		for i in range(len(cat)):
			this_cen = SkyCoord(cat['RA'].iloc[i]*u.deg,cat['Dec'].iloc[i]*u.deg)
			this_reg = CircleSkyRegion(center=this_cen,radius=cat['Radius_arcsec'].iloc[i]*u.arcsec)
			if label == True:
				this_reg.meta['label'] = str(i)
			this_reg.visual['fontsize'] = fontsize
			this_reg.visual['color'] = color
			regions.append(this_reg)

		return Regions(regions)

	def match_fovs(ref_cat,tar_cat,ref_wcs,tar_wcs):
		'''
		Find the minimally overlapping field of view between the reference and target catalogs.
		Remove sources outside this FoV.

		Keywords
		----------
		ref_cat : pandas DataFrame
			The "reference" catalog (e.g., MIRI), must have columns named 'RA' and 'Dec'
		tar_cat : pandas DataFrame
			The "target" catalog (e.g., others), must have columns named 'RA' and 'Dec'
		ref_wcs : astropy wcs object
			The wcs of the fits file used to contruct the reference catalog.
		tar_wcs : astropy wcs object
			The wcs of the fits file used to contruct the target catalog.

			
		Returns
		-------
		ref_cat : pandas DataFrame
			The updated reference catalog containing sources only within the overlapping FoVs
		tar_cat : pandas DataFrame
			The updated target catalog containing sources only within the overlapping FoVs
		ref_len : int
			The number of sources in the updated reference catalog
		tar_len : int
			The number of sources in the updated target catalog
		ndrop_ref_fov : int
			The number of sources removed from the original reference catalog that are outside the overlapping FoVs
		ndrop_tar_fov : int
			The number of sources removed from the original target catalog that are outside the overlapping FoVs
	 
		
		Notes
		-----
		Required packages: astropy, numpy, regions

		Author: Rebecca C. Levy (rlevy.astro@gmail.com)
		Last updated: 2025-07-29
		'''

		import numpy as np
		from astropy.wcs import WCS
		from astropy.coordinates import SkyCoord
		import astropy.units as u
		from regions import PolygonSkyRegion, Regions
		
		#convert to a SkyCoord object
		ref = SkyCoord(ref_cat['RA'],ref_cat['Dec'],unit=u.deg)
		tar = SkyCoord(tar_cat['RA'],tar_cat['Dec'],unit=u.deg)

		#determine smaller FoV
		ref_fov = PolygonSkyRegion(SkyCoord(ref_wcs.calc_footprint(),unit='deg'))
		tar_fov = PolygonSkyRegion(SkyCoord(tar_wcs.calc_footprint(),unit='deg'))
		min_fov = ref_fov.intersection(tar_fov) 

		#remove any target sources that are outside of the overlapping FoV
		ref_len_og = len(ref)
		tar_len_og = len(tar)

		#check whether each source is within the minimal FoV
		ref_in_min = min_fov.contains(ref,ref_wcs)
		tar_in_min = min_fov.contains(tar,ref_wcs)
		idx_ref_in_min = np.where(ref_in_min==True)[0]
		idx_tar_in_min = np.where(tar_in_min==True)[0]

		#keep only sources that are in the minimal FoV
		ref = ref[ref_in_min==True]
		tar = tar[tar_in_min==True]

		ref_cat_og = ref_cat.copy()
		tar_cat_og = tar_cat.copy()
		ref_cat = ref_cat.iloc[idx_ref_in_min]
		tar_cat = tar_cat.iloc[idx_tar_in_min]

		#calculate the new number of sources
		ref_len = len(ref_cat)
		tar_len = len(tar_cat)

		#calculate how many sources were removed from each catalog
		ndrop_ref_fov = ref_len_og - ref_len
		ndrop_tar_fov = tar_len_og - tar_len

		print('\t%i/%i reference sources within minimal FoV' %(ref_len,ref_len_og))
		print('\t%i/%i target sources within minimal FoV' %(tar_len,tar_len_og))

		return ref_cat, tar_cat, idx_ref_in_min, idx_tar_in_min, ndrop_ref_fov, ndrop_tar_fov

	def find_catalog_overlap(cat1,cat2,overlap_area_criterion=0):
		#match the catalogs based on the source sizes
		#cat1, cat2 should be a 3-tuple with (catalog,wcs,name)
		#overlap_area_criterion is the % overlapping area of the two sources to count as a match
		#the "reference catalog" is the catalog to match against; i.e., this is the catalog that will be searched against
		#the "target catalog" is the catalog to check for overlaps with ref_cat; 
		import numpy as np
		from astropy.wcs import WCS
		from astropy.coordinates import SkyCoord
		import astropy.units as u
		from regions import PolygonSkyRegion, Regions
		from JWSTTools import jwst_tools
		import sys

		#do a loop so each catalog can be the "reference" and the "target"
		cats = [[cat1,cat2],[cat2,cat1]]
		overlaps = []
		idx_overlaps = []
		fov_overlaps = []
		overlap_stats = []
		ref_cat_og = cat1
		tar_cat_og = cat2
		for j in range(len(cats)):
			ref_cat,ref_wcs,ref_name = cats[j][0]
			tar_cat,tar_wcs,tar_name = cats[j][1]

			print('------------------------------------')
			print('reference_catalog = %s' %ref_name)
			print('target_catalog = %s' %tar_name)
			print('------------------------------------')

			#convert to a SkyCoord object
			ref = SkyCoord(ref_cat['RA'],ref_cat['Dec'],unit=u.deg)
			tar = SkyCoord(tar_cat['RA'],tar_cat['Dec'],unit=u.deg)

			print('Removing sources outside of either FoV...')
			ref_cat, tar_cat, idx_ref_in_min, idx_tar_in_min, ndrop_ref_fov, ndrop_tar_fov \
				= jwst_tools.match_fovs(ref_cat,tar_cat,ref_wcs,tar_wcs)

			tar_len = len(tar_cat)
			ref_len = len(ref_cat)

			#loop over sources in the target catalog
			overlap_tar_cat = np.zeros((tar_len,),dtype=bool)
			idx_overlap_tar_cat = []
			for i in range(tar_len):
				sep = tar[i].separation(ref)
				#define overlapping sources as their centers being separated by less than the sum of their radii
				#i.e., they must overlap with one another within their respective radii
				if overlap_area_criterion == 0:
					criteria = sep.to(u.arcsec).value<=tar_cat['Radius_arcsec'].iloc[i]+ref_cat['Radius_arcsec']
				else:
					frac_overlap = jwst_tools.intersectionArea( 
									tar_cat['RA'].iloc[i]*u.deg,tar_cat['Dec'].iloc[i]*u.deg,tar_cat['Radius_arcsec'].iloc[i]*u.arcsec,
									ref_cat['RA'].values*u.deg,ref_cat['Dec'].values*u.deg,ref_cat['Radius_arcsec'].values*u.arcsec,
									ref_wcs)
					criteria = frac_overlap >= overlap_area_criterion

				if np.any(criteria):
					overlap_tar_cat[i]=True
					idx = np.where(criteria)[0][0]
					idx_overlap_tar_cat.append(idx)

			n_overlap = len(overlap_tar_cat[overlap_tar_cat==True])
			p_overlap = n_overlap/len(overlap_tar_cat)*100
			print('Source areas must overlap by %i%% to match' %(overlap_area_criterion*100))
			print('\t%1.f%% (%i/%i) of target sources have counterparts in reference'
				%(p_overlap,n_overlap,len(overlap_tar_cat)))

			overlap_stats.append([n_overlap,p_overlap/100,
								  tar_len,ndrop_tar_fov,tar_len+ndrop_tar_fov,
								  ref_len,ndrop_ref_fov,ref_len+ndrop_ref_fov])
			overlaps.append(overlap_tar_cat)
			idx_overlaps.append(idx_overlap_tar_cat)

			print('------------------------------------\n\n')

		cat1_in_cat2 = overlaps[0]
		cat2_in_cat1 = overlaps[1]
		idx_cat1_in_cat2 = idx_overlaps[0]
		idx_cat2_in_cat1 = idx_overlaps[1]
		cat1_fov = tar_cat.copy()
		cat2_fov = ref_cat.copy()

		return overlap_stats, [cat1_in_cat2,idx_cat1_in_cat2,cat1_fov], [cat2_in_cat1,idx_cat2_in_cat1,cat2_fov]
