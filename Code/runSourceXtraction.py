'''Extract point-like sources from image using sourceXtractor++
Based on runSourcExtraction.py by R.C.Levy used for M82 NIRCam star cluster catalog (Levy et al. 2024, ApJL)
Minor updates for SASP2025 for M82 MIRI F560W image
'''

import os
import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from regions import PixCoord,CirclePixelRegion,CircleSkyRegion,Regions
from JWSTTools import jwst_tools
import astropy.wcs as wcs
import warnings
import pandas as pd
from photutils.aperture import SkyCircularAperture,aperture_photometry
warnings.simplefilter('ignore', category=wcs.FITSFixedWarning)

def prune_overlapping_source(cat,wcs,psf_fwhm,frac_cut=0.75):
	#now drop any clusters which are signficantly contained within another cluster
	nc = len(cat[0])
	x2=cat[0]*u.deg
	y2=cat[1]*u.deg
	r2=psf_fwhm*u.arcsec*np.ones(x2.shape)
	for i in range(len(x2)):
		x1=x2[i]
		y1=y2[i]
		r1=r2[i]
		frac_overlaps = jwst_tools.intersectionArea(x1,y1,r1,x2,y2,r2,wcs)
		i_over = np.where(frac_overlaps >= frac_cut)[0]

		#replace the RA and Dec of overlapping sources with the average
		if i_over.size > 0:
			x2[i_over] = np.median(x2[i_over])
			y2[i_over] = np.median(y2[i_over])


	cat = np.array([x2,y2])
	cat = np.unique(cat,axis=1)
	print('\t%i clusters after overlap removal (-%i)' %(len(cat[0]),nc-len(cat[0])))
	nc = len(cat[0])
	return cat

def prune_overlapping_source_keep_brightest(cat,wcs,psf_fwhm,im,frac_cut=0.5):
	#now drop any clusters which are signficantly contained within another cluster
	nc = len(cat[0])
	x2=cat[0]*u.deg
	y2=cat[1]*u.deg
	r2=psf_fwhm*u.arcsec*np.ones(x2.shape)
	for i in range(len(x2)):
		x1=x2[i]
		y1=y2[i]
		r1=r2[i]
		frac_overlaps = jwst_tools.intersectionArea(x1,y1,r1,x2,y2,r2,wcs)
		i_over = np.where(frac_overlaps >= frac_cut)[0]

		#keep only the brightest clusters
		if i_over.size > 1:
			coords = SkyCoord(x2[i_over],y2[i_over])
			apers = SkyCircularAperture(coords,r1)
			apers_phot = [np.nansum(a.to_pixel(wcs=wcs).to_mask(method='center').multiply(im).ravel()) for a in apers]
			imax = np.argmax(apers_phot)
			x2[i_over] = x2[i_over[imax]]
			y2[i_over] = y2[i_over[imax]]


	cat = np.array([x2,y2])
	cat = np.unique(cat,axis=1)
	print('\t%i clusters after additional overlap removal (-%i)' %(len(cat[0]),nc-len(cat[0])))
	nc = len(cat[0])
	return cat

def prune_edge_source(cat,wcs,edge_buffer=10):
	nc = len(cat[0])
	#convert from wcs to pixels
	cat = wcs.all_world2pix(cat[0],cat[1],1,ra_dec_order=True)
	cat = np.array(cat).T

	naxis1 = wcs.array_shape[1]
	naxis2 = wcs.array_shape[0]

	cat_noedge = cat.copy()
	cat_noedge[cat[:,0]<edge_buffer,:] = np.nan
	cat_noedge[cat[:,1]<edge_buffer,:] = np.nan
	cat_noedge[cat[:,0]>naxis1-edge_buffer,:] = np.nan
	cat_noedge[cat[:,1]>naxis2-edge_buffer,:] = np.nan
	cat_noedge = cat_noedge[~np.isnan(cat_noedge)]
	cat_noedge = cat_noedge.reshape(int(len(cat_noedge)/2),2)

	#go back to wcs coords
	cat_noedge = wcs.all_pix2world(cat_noedge[:,0],cat_noedge[:,1],1,ra_dec_order=True)

	print('\t%i clusters after edge removal (-%i)' %(len(cat_noedge[0]),nc-len(cat_noedge[0])))
	nc = len(cat_noedge[0])

	return np.array(cat_noedge)

def cat2reg(cat,fits_header,radius):
	#cat should be in sky coordinates w/ units
	#radius needs unit

	wcs_frame = fits_header['RADESYS'].lower()
	
	regions = []
	for i in range(cat.shape[1]):
		this_cen = SkyCoord(cat[0][i],cat[1][i],frame=wcs_frame)
		this_reg = CircleSkyRegion(center=this_cen,radius=radius)
		this_reg.meta['label'] = str(i)
		this_reg.visual['fontsize'] = 5
		this_reg.visual['color'] = 'magenta'
		regions.append(this_reg)

	return Regions(regions)

def put_image_primary_hdu(fits_file_name,fits_file,out_path):
	#sourcextractor++ needs the image data to be in the Primary HDU of the fits file (extension 0)
	#but JWST data is in extention 1

	data = fits_file[1].data
	data[np.isnan(data)==True] = 0.
	header = fits_file[1].header
	new_fits_file = fits.HDUList([fits.PrimaryHDU(data=data,header=header)])
	new_fits_file.writeto(out_path+fits_file_name.replace('.fits','_imageonly.fits'), overwrite=True)

	return

def main(run_sourcextrator=True):
	out_path = '../Data/sourcextractor++/'

	#load the fits file and get the WCS info
	detection_filter = 'F560W'
	filter_info,_ = jwst_tools.MIRI_Filter_Info(this_filter=detection_filter)
	fits_file_name = jwst_tools.getFileNames(this_instrument='MIRI',this_filter=detection_filter)[0]
	fits_file = fits.open('../Data/JWST/'+fits_file_name)
	fits_header = fits_file[1].header
	det_im = fits_file[1].data
	wcs = WCS(fits_header)


	if run_sourcextrator == True:
		#run SourcExtractor++
		if os.environ['CONDA_DEFAULT_ENV'] != 'sourcex':
			print('Wrong environment! Please run\n\t$ conda activate sourcex\nand try again.')
			sys.exit()
		else:
			put_image_primary_hdu(fits_file_name,fits_file,out_path)
			os.system('sourcextractor++ --config-file ../Data/sourcextractor++/default.config')

	#open the catalog
	cat = np.loadtxt(out_path+'catalog.txt')
	cat_wcs = wcs.all_pix2world(cat[:,0],cat[:,1],1,ra_dec_order=True)

	#write catalog to region file in pixel and wcs units
	psf_fwhm_arcsec = filter_info['PSF_FWHM_arcsec']
	regions = cat2reg(cat_wcs*u.deg,fits_header,0.5*psf_fwhm_arcsec*u.arcsec)
	regions.write(out_path+'MIRI_F560Wsources_sourcextractor_regions_wcs.reg',format='ds9',overwrite=True)

	print('%i initial sources identified' %len(cat))


	#prune overlapping sources
	cat_pruned = prune_overlapping_source(cat_wcs,wcs,psf_fwhm_arcsec,frac_cut=0.75)

	#remove sources near the edges of the map
	cat_pruned = prune_edge_source(cat_pruned,wcs,edge_buffer=10)

	#prune slightly less overlapping sources
	cat_pruned = prune_overlapping_source_keep_brightest(cat_pruned,wcs,2*psf_fwhm_arcsec,det_im,frac_cut=0.5)


	#save cluster catalog
	np.savetxt(out_path+'catalog_pruned.txt',cat_pruned.T,fmt='%f',header='RA(deg) Dec(deg)')

	#save as region file
	regions_p = cat2reg(cat_pruned*u.deg,fits_header,0.5*psf_fwhm_arcsec*u.arcsec)
	regions_p.write(out_path+'MIRI_F560Wsources_sourcextractor_regions_wcs_pruned.reg',format='ds9',overwrite=True)

	#now save the catalog as a .csv with an ID number
	df = pd.DataFrame(cat_pruned.T,columns = ['RA','Dec'])
	df.to_csv(out_path+'MIRI_F560Wsources_sourcextractor.csv',index_label='ID')


	#copy to Catalogs folder
	os.system('cp '+out_path+'MIRI_F560Wsources_sourcextractor.csv ../Data/Catalogs/M82_ClusterCatalog_JWST-MIRI_sourcextractor.csv')


	return

if __name__ == '__main__':
	main()

