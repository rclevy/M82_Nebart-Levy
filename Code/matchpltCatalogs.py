this_script = __file__

import numpy as np
import pandas as pd
from JWSTTools import jwst_tools
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import ICRS,SkyCoord
import astropy.units as u
import astropy.wcs as wcs
from astropy.visualization import (ManualInterval, AsinhStretch, ImageNormalize, LogStretch, LinearStretch)
from regions import CircleSkyRegion, Regions
import warnings
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from reproject.mosaicking import find_optimal_celestial_wcs
from pltMIRIdata import add_scalebar
import sys
import datetime

warnings.simplefilter('ignore', category=wcs.FITSFixedWarning)
warnings.simplefilter('ignore', category=RuntimeWarning)

plt.rcParams['font.family']='serif'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'cm'


def get_catalog_info():
	#open the cluster catalogs
	path_to_catalogs = '../Data/Catalogs/'
	catalog_prefix = 'M82_ClusterCatalog_'
	catalog_names = ['JWST-MIRI_sourcextractor',
					'HST-ACS_Mayya2008',
					'HST-NICMOS_McCrady2003',
					'JWST-NIRCam_Levy2024',
					'SMA-850umcontinuum_JimenezDonaire']
	catalog_colors = ['xkcd:neon red',
					  'xkcd:azure',
					  'xkcd:cool green',
					  'xkcd:butter yellow',
					  'xkcd:pink']

	#open the fits files to get the wcs info
	fits_names = ['../Data/JWST/jw01701-o031_t007_miri_f560w-sub128_i2d.fits',
				'../Data/HST/hlsp_m82_hst_acs-wfc_all_f555w_v1_drz_rmSIP.fits',
				'../Data/HST/M82_HST_NICMOS_F190N_PaAlphaCont_mosaic.fits',
				'../Data/JWST/jw01701052001_nircam_sub640_f250m_jhat_i2d.fits',
				'../Data/SMA/m82_sma_cont.fits']
	return path_to_catalogs,catalog_prefix,catalog_names,catalog_colors,fits_names

def load_fits(fits_name):
	this_fits = fits.open(fits_name)
	if len(this_fits) > 1:
		ext=1
	else:
		ext=0
	this_fits = this_fits[ext]
	this_hdr = this_fits.header
	this_wcs = WCS(this_hdr)
	this_data = this_fits.data
	return this_fits, this_hdr, this_wcs, this_data

def load_catalog_from_reg(catalog_name):

	path_to_catalogs,catalog_prefix,_,_,_ = get_catalog_info()


	this_reg = Regions.read(path_to_catalogs+'region_files/'+catalog_prefix+catalog_name+'.reg')
	this_cat = jwst_tools.cat_from_reg(this_reg)

	if catalog_name == 'JWST-NIRCam_Levy2024':
		#use the F250M PSF HWHM instead of measured radius
		this_cat['Radius_arcsec'] = 2 * 0.085/2 + np.zeros_like(this_cat['Radius_arcsec'].values)
	elif catalog_name == 'JWST-MIRI_sourcextractor':
		this_cat['Radius_arcsec'] *= 2

	return this_reg, this_cat

def fmt_catalog_name_str(name,type='Obs+Inst'):
	if type == 'Obs+Inst':
		n = name.split('_')[0].replace('umcontinuum',r'$\mathrm{\mu m}$')
	elif type=='Obs':
		n = name.split('-')[0]
	elif type == 'Ref':
		n = name.split('_')[-1].replace('20','+20').replace('sourcextractor','Nebart+ip').replace('Donaire','Donaire+ip')
	elif type == 'Obs-or-Inst':
		if ('JWST' in name) or ('HST' in name):
			i = 1
		else:
			i = 0
		n = name.split('_')[0].split('-')[i]
	else:
		print('Invalid type specified.')
		n = name

	return n

def plot_image(ref_im,ref_wcs,ref_fits,norm_type='linear'):
	#plot the F560W image

	#we need to redefine as "reference" wcs so that north=up
	ref_wcs_nup,ref_shape_nup = find_optimal_celestial_wcs(ref_fits) 

	fig = plt.figure(2)
	plt.clf()
	ax = fig.add_subplot(projection=ref_wcs_nup)
	cmap='gray'
	if norm_type == 'asinh':
		norm = ImageNormalize(ref_im,
							 stretch=AsinhStretch(),
						 	interval=ManualInterval(vmin=0,vmax=3000))
	elif norm_type == 'log':
		norm = ImageNormalize(ref_im,
							  stretch=LogStretch(a=100),
							  interval=ManualInterval(vmin=0,vmax=5000))
	else:
		norm = ImageNormalize(ref_im,
							  stretch=LinearStretch(),
							  interval=ManualInterval(vmin=0,vmax=2000))
	ref_im[np.isnan(ref_im)==True]=0.
	im = ax.imshow(ref_im,origin='lower',cmap=cmap,norm=norm,transform=ax.get_transform(ref_wcs))
	ax.set_xlim(0,ref_shape_nup[1])
	ax.set_ylim(0,ref_shape_nup[0])
	mask_im = jwst_tools.mask_outside_footprint(ref_wcs_nup, ref_wcs, ref_shape_nup, ref_im.shape)
	ax.imshow(mask_im,origin='lower',cmap=cmap,transform=ax.get_transform(ref_wcs_nup))
	cb = plt.colorbar(im,label='Intensity (MJy/sr)',extend='max')
	ax.set_xlabel('R.A. (J2000)')
	ax.set_ylabel('Decl. (J2000)')
	ax.coords[0].set_major_formatter('hh:mm:ss') 
	ax.coords[0].set_separator((r'$^{\mathrm{h}}$',r'$^{\mathrm{m}}$',r'$^{\mathrm{s}}$'))
	ax.coords[0].set_ticklabel(exclude_overlapping=True)
	ax.coords[1].set_ticklabel(exclude_overlapping=True)
	ax.coords[0].display_minor_ticks(True) 
	ax.coords[1].display_minor_ticks(True)
	ax.coords[0].set_minor_frequency(4) 
	ax.coords[1].set_minor_frequency(3)
	_ = add_scalebar(ax,ref_wcs.to_header(),100,3.6E6,color='w')
	return fig, ax, ref_wcs_nup, norm_type

def plot_catalogs(fig,ax,ref_wcs_nup,catalog_regions,catalog_colors,catalog_names,idx_overlaps):
	#now add all the catalogs
	leg = []
	n_reg = len(catalog_regions)

	alpha=1.0
	lw=0.25

	for i in range(n_reg):
		this_color = catalog_colors[i]
		if idx_overlaps == None:
			this_overlap = None
		else:
			this_overlap = idx_overlaps[i]
		
		if i == 0:
			zorder=10
		else:
			zorder=8
		for j in range(len(catalog_regions[i])):
			if np.all(this_overlap != None):
				if j in this_overlap:
					alpha=1.0
					lw=0.5
				else:
					lw=0.25
					# alpha=0.33
					alpha=0.5
			this_reg = catalog_regions[i][j].to_pixel(ref_wcs_nup)
			this_reg.plot(ec=this_color,alpha=alpha,lw=lw)#,zorder=zorder)
		leg_str = catalog_names[i].split('_')[0].replace('-','/').replace('umcontinuum',r'$\mathrm{\mu m}$')
		# leg_str = catalog_names[i].replace('_',' ').replace('-',' ').replace('umcontinuum',r'$\mathrm{\mu m}$')
		l,=ax.plot(np.nan,np.nan,'o',color=this_color,mfc='None',label=leg_str)

		leg.append(l)

	if len(leg) > 3:
		leg = [leg[i] for i in [1,2,3,0,4]]
	plt.legend(leg,[ll.get_label() for ll in leg],
		loc='upper right',fontsize=plt.rcParams['font.size']-2,
		labelcolor='w',facecolor='k',edgecolor='w')
	return fig,ax

def make_uber_catalog(oac=0.5):
	#combine all of the catalogs together to make a giant cluster catalog
	#denote the original source of each cluster and whether it overlaps with others
	path_to_catalogs,catalog_prefix,catalog_names,catalog_colors,fits_names = get_catalog_info()
	catalog_suffix = 'AreaOverlap%i' %(oac*100)

	#sort by wavelength
	catalog_names = [catalog_names[i] for i in [1,2,3,0,4]]
	fits_names = [fits_names[i] for i in [1,2,3,0,4]]

	#format catalog references
	cat_refs = [c.split('_')[-1].replace('20','+20').replace('sourcextractor','Nebart+ip').replace('Donaire','Donaire+ip') for c in catalog_names]

	#set column order
	fmt_cat_refs = ['in_'+c for c in cat_refs]
	col_order = ['ID','RA','Dec','Radius_arcsec','Reference']+fmt_cat_refs

	uber_cat = []

	for i in range(len(catalog_names)):
		cat_name = catalog_names[i]

		#open the catalog csv
		cat = pd.read_csv(path_to_catalogs+catalog_suffix+'/'+catalog_prefix+cat_name+'_CrossMatches_'+catalog_suffix+'.csv')

		#add the reference for this catalog
		cat['Reference'] = cat_refs[i]

		#rename the columns
		for j in range(len(fmt_cat_refs)):
			col_to_replace = catalog_names[j].split('_')[1]
			cat.rename(columns={col_to_replace:fmt_cat_refs[j]},inplace=True)

		#set overlaps with itself to True
		cat[fmt_cat_refs[i]] = True

		#reorder columns
		cat = cat[col_order]

		uber_cat.append(cat)

	#now join into one long dataframe
	df = pd.concat(uber_cat)
	df.to_csv(path_to_catalogs+catalog_suffix+'/'+catalog_prefix.replace('_C','_AllC').replace('g_','gs_')+'CrossMatches_'+catalog_suffix+'.csv',index=False)

	#write a txt file that describes this catalog
	with open(path_to_catalogs+catalog_suffix+'/'+catalog_prefix.replace('_C','_AllC').replace('g_','gs_')+'CrossMatches_'+catalog_suffix+'.txt','w') as f:
		f.write(catalog_prefix.replace('_C','_AllC').replace('g_','gs_')+'CrossMatches_'+catalog_suffix+'.txt\n')
		f.write('Written on %s\n\n' %datetime.datetime.now())
		f.write('The combined star cluster candidate catalogs for M82. The cross-matching assumes that cluster areas must overlap by >%i%%. Cross-matches are only reported within the MIRI subarray FoV, the smallest among the cluster catalogs.\n\n' %(oac*100))
		f.write('Data Sources:\n')
		f.write('\tMayya+2008 \t HST-ACS\n\tMcCrady+2003 \t HST-NICMOS\n\tLevy+2024 \t JWST-NIRCam\n\tNebart+in \t JWST-MIRI\n\tJimenezDonaire+ip \t SMA-850micron\n\n')
		f.write('Column Descriptions:\n')
		f.write('\tID \t The cluster identification in the original catalog\n')
		f.write('\tRA \t The J2000 right ascension of the cluster center in the original catalog. For HST-based catalogs, the WCS coordinates have been corrected.\n')
		f.write('\tDec \t The J2000 declination of the cluster center in the original catalog. For HST-based catalogs, the WCS coordinates have been corrected.\n')
		f.write('\tRadius_arcsec \t The cluster radius (in arcsec) used for the cross-matching.\n')
		f.write('\tReference \t The original reference for this star cluster candidate.\n')
		f.write('\tin_Mayya+2008 \t Whether the given cluster overlaps with source(s) in the Mayya+2008 catalog. Values are True (overlap found), False (no overlap found), or blank (cluster outside MIRI FoV). Matches with the Reference are always True.\n')
		f.write('\tin_McCrady+2003 \t The same as above, but for the McCrady+2003 catalog.\n')
		f.write('\tin_Levy+2024 \t The same as above, but for the Levy+2024 catalog.\n')
		f.write('\tin_Nebart+ip \t The same as above, but for the Nebart+ip catalog.\n')
		f.write('\tin_JimenezDonaire+ip \t The same as above, but for the JimenezDonaire+ip catalog.\n')

	return df

def plot_supervenn(oac=0.5):
	from supervenn import supervenn, make_sets_from_chunk_sizes

	catalog_suffix = 'AreaOverlap%i' %(oac*100)
	uber_cat = make_uber_catalog(oac=oac)
	_,_,catalog_names,catalog_colors,_ = get_catalog_info()
	uc = uber_cat.dropna(axis=0).reset_index()
	uc = uc.iloc[:,[-5,-4,-3,-2,-1]]
	cn = [b.split('_')[0].replace('-','/').replace('umcontinuum','$\\mu$m') for b in catalog_names]
	cn = [cn[i] for i in [1,2,3,0,4]]
	cols = [catalog_colors[i] for i in [1,2,3,0,4]]
	uc = uc.rename(columns={old:new for (old,new) in zip(uc.keys().values.tolist(),cn)})
	uc['size'] = np.ones((len(uc),)).astype(int)

	sets,labels = make_sets_from_chunk_sizes(uc)
	sets.reverse()
	labels.reverse()
	cols.reverse()

	plt.figure(1,figsize=(20,6))
	plt.clf()
	ax = plt.gca()
	supervenn(sets,labels,ax=ax,
		widths_minmax_ratio=0.05,min_width_for_annotation=-1,
		color_cycle=cols,bar_alpha=1,
		rotate_col_annotations=False,alternating_background=False,
		side_plots=False,col_annotations_area_height=0.4,
		chunks_ordering='occurrence',
		fontsize=plt.rcParams['font.size']+2)#,reverse_chunks_order=False)
	plt.ylabel('Catalogs')
	plt.xlabel('Number of YMC Candidates in each Catalog Combination')
	plt.savefig('../Plots/M82_AllClusterCatalogs_Supervenn_'+catalog_suffix+'.pdf',
		bbox_inches='tight',metadata={'Creator':this_script})
	plt.close()

	return

def plot_histogram(oac=0.5):
	import matplotlib as mpl
	mpl.rcParams['hatch.linewidth'] = 8.0 

	path_to_catalogs,catalog_prefix,catalog_names,catalog_colors,_ = get_catalog_info()


	#open all the overlap files
	def load_olaps(oac):
		df = []
		catalog_suffix = 'AreaOverlap%i' %(oac*100)
		pre = path_to_catalogs+catalog_suffix+'/'+catalog_prefix.replace('g_','gs_')+'OverlapStats_'
		suff = '_'+catalog_suffix+'.csv'
		fnames = [pre+s+suff for s in catalog_names]
		[df.append(pd.read_csv(f)) for f in fnames]
		df = pd.concat(df)
		df.drop_duplicates(inplace=True)
		df.drop(df.columns[[3,5,6,7,8,9]],axis=1,inplace=True)
		df = df.sort_values(['Target_Catalog','Reference_Catalog'],axis=0)
		catalog_names.sort()

		#get overlaps and unique numbers
		n_a = df.drop_duplicates(subset='Target_Catalog')['N_tar-fov'].values
		n_a_in_b = [df[df['Target_Catalog']==c]['N_tar_in_ref'].values for c in catalog_names]
		n_b_in_a = [df[df['Reference_Catalog']==c]['N_tar_in_ref'].values for c in catalog_names]
		return n_a, n_a_in_b, n_b_in_a, catalog_names, df, catalog_suffix

	n_a, n_a_in_b, n_b_in_a, catalog_names, df, catalog_suffix = load_olaps(oac)

	overlap_ave = np.mean(np.array([n_a_in_b,n_b_in_a]),axis=0)
	overlap_min = np.min(np.array([n_a_in_b,n_b_in_a]),axis=0)
	overlap_max = np.max(np.array([n_a_in_b,n_b_in_a]),axis=0)

	n_a_only = np.abs(n_a - np.sum(overlap_ave,axis=1))
	n_a_only_min = np.abs(n_a - np.sum(overlap_max,axis=1))
	n_a_only_max = np.abs(n_a - np.sum(overlap_min,axis=1))

	#set up the histogram
	order = [4,2,3,1,0]
	catalog_colors = [catalog_colors[i] for i in [4,0,3,2,1]]
	bin_names_only = [catalog_names[i] for i in order]
	n_only = [n_a_only[i] for i in order]
	n_only_min = [n_a_only_min[i] for i in order]
	n_only_max = [n_a_only_max[i] for i in order]
	bin_names_overlap = bin_names_only[1:]

	bn = [b.split('_')[0].replace('-','/').replace('continuum','') for b in bin_names_only]
	bn.reverse()
	set_nums = np.concatenate((n_a_only.reshape(len(n_a_only),1),overlap_ave),axis=1)

	tmp = []
	for i in range(len(bn)):
		tmp.append([[b+str(i) for (b,i) in zip(bn[i],range(int(set_nums[0,i])))]])

	n_unique_total,_ = count_unique_sources(oac=oac)


	x = []
	y = []
	eyu = []
	eyl = []
	c = []
	h = []
	hc = []

	for i in range(len(bin_names_overlap)):

		n1 = df[(df['Target_Catalog']==bin_names_only[i]) & (df['Reference_Catalog']==bin_names_overlap[i])]['N_tar_in_ref'].values
		n2 = df[(df['Reference_Catalog']==bin_names_only[i]) & (df['Target_Catalog']==bin_names_overlap[i])]['N_tar_in_ref'].values
		n_overlap = np.mean([n1,n2])
		n_overlap_min = np.min([n1,n2])
		n_overlap_max = np.max([n1,n2])

		# n10 = df0[(df0['Target_Catalog']==bin_names_only[i]) & (df0['Reference_Catalog']==bin_names_overlap[i])]['N_tar_in_ref'].values
		# n20 = df0[(df0['Reference_Catalog']==bin_names_only[i]) & (df0['Target_Catalog']==bin_names_overlap[i])]['N_tar_in_ref'].values
		# n_overlap_min0 = np.min([n10,n20])
		# n_overlap_max0 = np.max([n10,n20])

		x0 = fmt_catalog_name_str(bin_names_only[i],type='Obs-or-Inst')
		x1 = fmt_catalog_name_str(bin_names_overlap[i],type='Obs-or-Inst')
		x.append(x0+'\nonly')
		x.append(x0+'\n&\n'+x1)
		y.append(n_only[i])
		eyu.append(n_only_max[i])
		eyl.append(n_only_min[i])
		y.append(n_overlap)
		eyu.append(n_overlap_max)
		eyl.append(n_overlap_min)
		c.append(catalog_colors[i])
		c.append(catalog_colors[i])
		h.append(None)
		h.append('\\')
		hc.append('None')
		hc.append(catalog_colors[i+1])

		if i == len(bin_names_overlap)-1:
			x0 = fmt_catalog_name_str(bin_names_only[i+1],type='Obs-or-Inst')
			x.append(x0+'\nonly')
			y.append(n_only[i+1])
			eyu.append(n_only_max[i+1])
			eyl.append(n_only_min[i+1])
			c.append(catalog_colors[i+1])
			h.append(None)
			hc.append('None')


	y = np.array(y)
	eyu = np.array(eyu)
	eyl = np.array(eyl)

	eyu = np.abs(eyu-y)
	eyl = np.abs(eyl-y)

	fig = plt.figure(1,figsize=(9,3))
	plt.clf()
	plt.bar(x,y,
		log=True,
		ec=hc,
		fc=c,
		hatch=h,
		zorder=9)
	#add error bars
	b=plt.bar(x,y,yerr=(eyl,eyu),
		log=True,
		ec='k',
		lw=1.5,
		fc='None',
		zorder=10,
		capsize=3)

	b_labels = ['%.1f%%' %(yy/n_unique_total*100) for yy in y]
	plt.bar_label(b,labels=b_labels,fontsize=plt.rcParams['font.size']-4)

	plt.ylim(bottom=1,top=2000)

	plt.axhline(n_unique_total,color='k',lw=2)
	plt.text(0.005,0.92,'N$_{\\mathrm{unique}}$=%i' %(n_unique_total),
		fontsize=plt.rcParams['font.size']-3,ha='left',va='bottom',
		transform=plt.gca().transAxes)

	plt.xticks(fontsize=plt.rcParams['font.size']-2,)
	plt.gca().tick_params(axis='y',which='both',direction='in',right=True)
	plt.grid(alpha=0.5,axis='y')
	plt.grid(alpha=0.1,axis='y',which='minor')

	plt.ylabel('Number of massive\nstar cluster candidates')
	plt.xlabel('Relative age $\\rightarrow$',labelpad=10)
	for (x,t,ha) in zip([0,1],['(younger)','(older)'],['left','right']):
		plt.text(x,-0.2,t,
			fontsize=plt.rcParams['font.size']-3,
			ha=ha,va='top',transform=plt.gca().transAxes)

	plt.savefig('../Plots/M82_AllClusterCatalogs_Histogram_'+catalog_suffix+'.pdf',
		bbox_inches='tight',metadata={'Creator':this_script})
	plt.close()

	return

def count_unique_sources(oac=0.5):
	#count unique clusters using an extension of set theory
	path_to_catalogs,catalog_prefix,catalog_names,catalog_colors,fits_names = get_catalog_info()
	catalog_suffix = 'AreaOverlap%i' %(oac*100)

	#open all the overlap files
	df = []
	pre = path_to_catalogs+catalog_suffix+'/'+catalog_prefix.replace('g_','gs_')+'OverlapStats_'
	suff = '_'+catalog_suffix+'.csv'
	fnames = [pre+s+suff for s in catalog_names]
	[df.append(pd.read_csv(f)) for f in fnames]
	df = pd.concat(df)
	df.drop_duplicates(inplace=True)
	df.drop(df.columns[[3,5,6,7,8,9]],axis=1,inplace=True)
	n_a_in_b = df['N_tar_in_ref'].iloc[::2].values
	n_b_in_a = df['N_tar_in_ref'].iloc[1::2].values
	n_a = df.drop_duplicates(subset='Target_Catalog')['N_tar-fov'].values
	print('%i entries from all catalogs within FoV' %np.sum(n_a))

	overlap_max = np.max([n_a_in_b,n_b_in_a],axis=0)
	overlap_min = np.min([n_a_in_b,n_b_in_a],axis=0)
	overlap_ave = np.mean([n_a_in_b,n_b_in_a],axis=0)

	n_unique_ave = np.sum(n_a) - np.sum(overlap_ave)
	n_unique_min = np.sum(n_a) - np.sum(overlap_max)
	n_unique_max = np.sum(n_a) - np.sum(overlap_min)

	n_unique_total = n_unique_ave.copy()
	en_unique_total = n_unique_ave - n_unique_min


	print('\n\t##########################################################\n')
	print('\tM82 has %i +/- %i unique massive star cluster candidates\n\twithin the MIRI FoV and assuming %i%% overlap\n'
		%(n_unique_total,en_unique_total,oac*100))
	print('\t##########################################################\n')


	return n_unique_total, en_unique_total

def plot_all_catalog_overlaps(oac,uber_cat):
	cats = np.unique(uber_cat['Reference'].values)
	_,_,catalog_names,catalog_colors,fits_names = get_catalog_info()
	catalog_names = [catalog_names[i] for i in [1,2,3,0,4]]
	catalog_colors = [catalog_colors[i] for i in [1,2,3,0,4]]
	cats = [cats[i] for i in [2,3,1,4,0]]

	ref_fits, ref_hdr, ref_wcs, ref_im = load_fits(fits_names[0])
	fig,ax,ref_wcs_nup,norm_type = plot_image(ref_im,ref_wcs,ref_fits)


	for j in range(len(cats)):
		sub_uc = uber_cat[uber_cat['Reference']==cats[j]]
		sub_uc = sub_uc.drop(columns='in_'+cats[j])
		cms = sub_uc.iloc[:,[-4,-3,-2,-1]]
		idx_overlap = []
		for i in range(len(sub_uc)):
			if np.any(cms.iloc[i,:] == True):
				idx_overlap.append(i)
		reg = jwst_tools.reg_from_cat(sub_uc)
		fig,ax = plot_catalogs(fig,ax,ref_wcs_nup,[reg],[catalog_colors[j]],[catalog_names[j]],[idx_overlap])


	leg = ax.legend(loc='upper right',fontsize=plt.rcParams['font.size']-2,
	labelcolor='w',facecolor='k',edgecolor='w')
	
	catalog_suffix = 'AreaOverlap%i' %(oac*100)		
	plt.savefig('../Plots/M82_MIRI-F560W_AllClusterCatalogs_CrossMatches_'+catalog_suffix+'.pdf',bbox_inches='tight',metadata = {'Creator':this_script})
	plt.close('all')
	return

def main(run_only_post = False):

	path_to_catalogs,catalog_prefix,catalog_names,catalog_colors,fits_names = get_catalog_info()


	#define the criteria for "overlapping"
	#based on the overlapping source areas
	#default = 0 which means just barely touching
	# overlap_area_criterion = [0,0.5]
	overlap_area_criterion = [0.5]


	#always use the MIRI FoV as the reference since it's the smallest
	_,_,fov_ref_wcs,_ = load_fits(fits_names[0])


	if run_only_post == False:
		for ii in range(len(catalog_names)):

			ref_name = catalog_names[ii]
			ref_reg, ref_cat = load_catalog_from_reg(ref_name)
			ref_fits, ref_hdr, ref_wcs, ref_im = load_fits(fits_names[ii])

			#print('%s: %i' %(ref_name,len(ref_reg)))

			if ii == 0:
				do_plot = True
			else:
				do_plot = False

			print_to_terminal = False

			#loop over different "overlap" criteria
			for oac in overlap_area_criterion:
				catalog_suffix = 'AreaOverlap%i' %(oac*100)

				if print_to_terminal == False:
					f = open(path_to_catalogs+catalog_suffix+'/logs/'+this_script.split('/')[-1].replace('.py','')+'-output_'+ref_name+'.txt','w')
					sys.stdout = f
					print('Output on:\n%s\n\n' %datetime.datetime.now())

				#set up arrys to hold info
				ref_in_target = []
				# target_in_ref = []
				idx_ref_in_target = []
				idx_target_in_ref = []
				# ref_coords = []
				catalog_regions = [ref_reg]
				key_list = ['Target_Catalog','Reference_Catalog',
							'N_tar_in_ref','Frac_tar_in_ref',
							'N_tar-fov','N-drop_tar','N_tar-og',
							'N_ref-fov','N-drop_ref','N_ref-og']
				dtype_list = [str,str,np.float64,np.float64,np.float64,np.float64,np.float64,np.float64,np.float64,np.float64]
				dict_overlap_stats = {k: np.array([]).astype(d) for k,d in zip(key_list,dtype_list)}

				#now loop over the other "comparison catalogs"
				for jj in range(len(catalog_names)):
					target_name = catalog_names[jj]
					if target_name != ref_name:
						
						target_reg, target_cat = load_catalog_from_reg(target_name)
						target_fits, target_hdr, target_wcs, target_im = load_fits(fits_names[jj])

						#match the catalogs
						#this removes sources outside the ref_name FoV
						overlap_stats, overlap_info_ref_cat, overlap_info_target_cat = \
							jwst_tools.find_catalog_overlap((ref_cat,fov_ref_wcs,ref_name),
															(target_cat,fov_ref_wcs,target_name),
															overlap_area_criterion=oac)

						ref_target_names = [[target_name,ref_name],[ref_name,target_name]]
						for k,j in zip(key_list[0:2],ref_target_names):
							dict_overlap_stats[k] = np.concatenate((dict_overlap_stats[k],j))
						for k,j in zip(key_list[2:],np.array(overlap_stats).T):
							dict_overlap_stats[k] = np.concatenate((dict_overlap_stats[k],j))


						ref_coords_overlap = overlap_info_ref_cat[2].reset_index(drop=True)
						ref_in_target.append(overlap_info_target_cat[0])
						# ref_in_target.append(overlap_info_ref_cat[0])
						idx_ref_in_target.append(overlap_info_ref_cat[1])
						# target_in_ref.append(overlap_info_target_cat[0])
						idx_target_in_ref.append(overlap_info_target_cat[1])

						target_reg_fov = jwst_tools.reg_from_cat(overlap_info_target_cat[2])
						catalog_regions.append(target_reg_fov)

						# # #write out the overlap with the aux catalog
						# ref_coords.append(overlap_info_ref_cat[2])
						# # dff[ref_name] = overlap_info_ref_cat[0]
						# # dff.to_csv(path_to_catalogs+catalog_suffix+'/'+catalog_prefix+target_name+'_CrossMatch_'+catalog_suffix+'.csv',index=False)

				#write the MIRI catalog and cross-matches to a file
				cn = catalog_names.copy()
				cn.remove(ref_name)
				target_names_short = [t.split('_')[-1] for t in cn]
				df_overlaps = pd.DataFrame(data=np.array(ref_in_target).T,columns=target_names_short)
				df_overlaps = pd.concat([ref_coords_overlap,df_overlaps],axis=1)
				#combine back with full reference catalog

				df = ref_cat.merge(df_overlaps,left_on='ID',right_on='ID',how='left')
				df = pd.concat([ref_cat,df],axis=1)
				df.to_csv(path_to_catalogs+catalog_suffix+'/'+catalog_prefix+ref_name+'_CrossMatches_'+catalog_suffix+'.csv',index=False)

				#write the overlap stats to a file
				df_stats = pd.DataFrame.from_dict(dict_overlap_stats)
				df_stats.to_csv(path_to_catalogs+catalog_suffix+'/'+catalog_prefix.replace('g_','gs_')+'OverlapStats_'+ref_name+'_'+catalog_suffix+'.csv',index=False)

				if print_to_terminal == False:
					f.close()
					sys.stdout = sys.__stdout__

				if do_plot == True:
					# PLOT
					miri_matches = np.unique(np.concatenate(idx_ref_in_target).ravel())
					idx_overlaps = idx_target_in_ref
					idx_overlaps.insert(0,miri_matches)

					idx_miri_overlaps = idx_ref_in_target
					idx_miri_overlaps.insert(0,miri_matches)


					#plot the F560W image
					fig,ax,fov_ref_wcs_nup,norm_type = plot_image(ref_im,fov_ref_wcs,ref_fits)
					plt.savefig('../Plots/M82_MIRI-F560W_'+norm_type+'.pdf',bbox_inches='tight',metadata = {'Creator':this_script})

					# #add all the cluster catalogs, showing overlaps with MIRI
					# fig,ax = plot_catalogs(fig,ax,fov_ref_wcs_nup,
					# 					   catalog_regions,
					# 					   catalog_colors,
					# 					   catalog_names,
					# 					   idx_overlaps)
					# plt.savefig('../Plots/M82_MIRI-F560W_AllClusterCatalogs_CrossMatches_'+catalog_suffix+'.pdf',bbox_inches='tight',metadata = {'Creator':this_script})
					# plt.close()

					# #plot each indivdual catalog and overlaps with MIRI
					# for jj in range(len(catalog_names)):
					# 	target_name = catalog_names[jj]
					# 	if target_name != ref_name:
					# 		fig,ax,fov_ref_wcs_nup,norm_type = plot_image(ref_im,fov_ref_wcs,ref_fits)
					# 		fig,ax = plot_catalogs(fig,ax,fov_ref_wcs_nup,
					# 							  [catalog_regions[0],catalog_regions[jj]],
					# 							  [catalog_colors[0],catalog_colors[jj]],
					# 							  [catalog_names[0],catalog_names[jj]],
					# 							  [idx_miri_overlaps[jj],idx_overlaps[jj]])
					# 		plt.savefig('../Plots/M82_MIRI-F560W_'+ref_name+'_'+target_name+'_CrossMatches_'+catalog_suffix+'.pdf',bbox_inches='tight',metadata = {'Creator':this_script})
					# 		plt.close()

					#plot just the MIRI catalog
					fig,ax,fov_ref_wcs_nup,norm_type = plot_image(ref_im,fov_ref_wcs,ref_fits)
					fig,ax = plot_catalogs(fig,ax,fov_ref_wcs_nup,[catalog_regions[0]],[catalog_colors[0]],[catalog_names[0]],None)
					plt.savefig('../Plots/M82_MIRI-F560W_ClusterCatalog_'+ref_name+'.pdf',bbox_inches='tight',metadata = {'Creator':this_script})
					plt.close()



	sys.stdout = sys.__stdout__
	for oac in overlap_area_criterion:
		# n_unique_total, en_unique_total = count_unique_sources(oac=oac)
		plot_histogram(oac=oac)
		plot_supervenn(oac=oac)
		uber_cat = make_uber_catalog(oac=oac)
		plot_all_catalog_overlaps(oac,uber_cat)
		

	return 



if __name__ == '__main__':
	main(run_only_post=True)


