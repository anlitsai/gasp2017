#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 05:22:35 2019

@author: altsai
"""

import os
import sys
import shutil
import numpy as np
import csv
import time
import math
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord  # High-level coordinates
#from astropy.coordinates import ICRS, Galactic, FK4, FK5 # Low-level frames
#from astropy.coordinates import Angle, Latitude, Longitude  # Angles
#from astropy.coordinates import match_coordinates_sky
from astropy.table import Table
from photutils import CircularAperture
from photutils import SkyCircularAperture
from photutils import aperture_photometry
from photutils import CircularAnnulus
from photutils import SkyCircularAnnulus
# https://photutils.readthedocs.io/en/stable/aperture.html
#from phot import aperphot
# http://www.mit.edu/~iancross/python/phot.html

import matplotlib.pyplot as plt
import matplotlib.axes as ax
from astropy.io import fits
from astropy.wcs import WCS
#from photutils import DAOStarFinder
#from astropy.stats import mad_std
# https://photutils.readthedocs.io/en/stable/getting_started.html

from numpy.polynomial.polynomial import polyfit
from astropy.stats import sigma_clipped_stats
from photutils.psf import IterativelySubtractedPSFPhotometry
from statistics import mode
from astropy.visualization import simple_norm
from photutils.utils import calc_total_error
from astropy.stats import mad_std
import matplotlib.gridspec as gridspec

import julian
import datetime


'''
# 3C345-20190714@135841-R_Astrodon_2018_calib.fits
#fits_root=input('Please input the file name of the fitsinstmage: ').split('.',-1)[0].split('_calib',-1)[0]
fits_root='3C345-20190822@130653-R_Astrodon_2018'
fits_calib=fits_root+'_calib.fits'
fits_ori=fits_root+'.fts'
#print(fits_root)
#print(fits_calib)
#print(fits_ori)

date=fits_root.split('@',-1)[0].split('-',-1)[-1]
print(date)
print(type(date))
yearmonth=date[0:6]
'''
#file_info='gasp_target_fitsheader_info.txt'
file_info='gasp_target_fitsheader_info_exclude_baddata_2017.txt'
print('... will read '+file_info+' ...')
df_info=pd.read_csv(file_info,delimiter='|')
print(df_info)

#obj_name='3C66A'
#obj_name='3C345'
#obj_name='3C279'
#obj_name='3C371'
#obj_name='3C454-3'
#obj_name='ES2344+514'
#obj_name='L-Lacertae'
#obj_name='4C38-41'
#obj_name='4C71-07'
#obj_name='4C29-45'
#obj_name='Mkn421'
obj_name='Mkn501'
#obj_name='AO0235+16'

#obj_name=sys.argv[1]


dir_obj='Rmag_InstMag/annu_w1/'+obj_name+'/'
print('... will generate files in: ./'+dir_obj)
if os.path.exists(dir_obj):
    shutil.rmtree(dir_obj)
os.makedirs(dir_obj,exist_ok=True)

dir_refstar='RefStar/'
#file_refstar='gasp_refStar_radec_annu_test.txt'
file_refstar='gasp_refStar_radec_annu.txt'

df_refstar=pd.read_csv(file_refstar,sep='|')
#print(df_refstar)
idx_refstar=df_refstar[df_refstar['ObjectName']==obj_name].index.tolist()
#print(idx_refstar)
n_refstar=len(idx_refstar)
refstarID=['']*n_refstar
#print('number of reference stars : ',n_refstar)
ra_deg=np.array([0.]*n_refstar)
dec_deg=np.array([0.]*n_refstar)
#ra_pix=np.array([0.]*n_refstar)
#dec_pix=np.array([0.]*n_refstar)
radec_deg=np.array([[0.,0.]]*n_refstar)
#radec_pix=np.array([[0.,0.]]*n_refstar)
ra_hhmmss=['']*n_refstar
dec_ddmmss=['']*n_refstar

rmag=np.array([0.]*n_refstar)
rmag_err=np.array([0.]*n_refstar)
rmag_uplims=np.array([0.]*n_refstar,dtype=bool)
#mag_instrument=np.array([0.]*n_refstar)
#mag_err=np.array([0.]*n_refstar)
#snr=np.array([0.]*n_refstar)
r_circle=np.array([0.]*n_refstar)
r_inner=np.array([0.]*n_refstar)
r_outer=np.array([0.]*n_refstar)


#positions=['']*n_refstar
aperture=['']*n_refstar
aper_annu=['']*n_refstar
aperture_pix=['']*n_refstar
aper_annu_pix=['']*n_refstar


font_small=8
font_middle=10
font_large=12
colors=['red','orange','green','blue','purple','brown','olive','pink']
markers=['o','v','^','s','D','<','>','p']

print('-----------------------')

# here is for test
#ID_selected=np.array([1276,1277,1371,1372])   # 3C345
#ID_selected=np.array([810,811,1334,1335])   # ES2344+514
#idx_fitsheader=np.array([809,810,1333,1334,1575,1576])   # ES2344+514
#ID_selected=np.array([814,815,1336,1337])  # L-Lacertae
#ID_selected=np.array([1330,1331,1476,1477,1478,1479])   # 3C66A
#ID_selected=np.array([810,811,1040,1041])   # ES2344+514
#ID_selected=np.array([774,775,776,777])   # 4C38-41
#ID_selected=np.array([1036,1037])   # ES2344+514
#ID_selected=np.array([1037,1038])   # ES2344+514
#ID_selected=np.array([1436,1437,1438,1439])   # ES2344+514
#ID_selected=np.array([988,989]) # 4C29-45
#ID_selected=np.array([963,964,991,992])   # 3C279
#ID_selected=np.array([687,688,995,996])   # 4C71-07
#ID_selected=np.array([689,690,997,998])   # Mkn421
#ID_selected=np.array([1527,1528,1533,1534])  #AO0235+16

#idx_fitsheader=ID_selected-1
idx_fitsheader=df_info[df_info['Object']==obj_name].index
print(idx_fitsheader)
#obj_name=df_info['Object'][idx_fitsheader]

fits_ori=df_info['Filename'][idx_fitsheader]
print(fits_ori)

#sys.exit(0)

#=======================

df_baddata_note=pd.read_csv('bad_data_note_2017.txt',sep='\t')
#print(df_baddata_note)

df_cannotfit=df_baddata_note.loc[(df_baddata_note['NoteIdx']==3)].reset_index(drop=True)
n_cannotfit=len(df_cannotfit)
print('... there are',n_cannotfit,'fits can not be fitted ...')
list_cannotfit=df_cannotfit['filename'].tolist()
print(list_cannotfit)


#idx_fitsheader
#del idx_fitsheader[idx_cannotfit]

n_idx=len(idx_fitsheader)

idx_cannotfit=[]
print('before drop bad fits: idx_fitsheader =',idx_fitsheader)
for i in range(n_idx):
    j=idx_fitsheader[i]
    fitsname=fits_ori[j]
#    print(j,fitsname)
    if fitsname in list_cannotfit:
        idx_cannotfit.append(i)
#        print(idx_cannotfit)
        #ignore=ignore+1
        print(i,j,fitsname,'... bad fits')
#        print('skip',j,fitsname,'... bad fits ...')
        print()
    else:
        print(i,j,fitsname)


idx_fitsheader_canfit=np.delete(idx_fitsheader,idx_cannotfit)
print(idx_fitsheader_canfit)
ID_canfit=idx_fitsheader_canfit+1
print(ID_canfit)
idx_fitsheader_cannotfit=idx_fitsheader[idx_cannotfit]
ID_cannotfit=idx_fitsheader_cannotfit+1
print(ID_cannotfit)
n_idx_canfit=len(idx_fitsheader_canfit)
n_cannotfit=len(idx_cannotfit)
print()
#sys.exit(0)     

n_idx=n_idx_canfit
idx_fitsheader=idx_fitsheader_canfit
#=======================

#n_idx=len(idx_fitsheader)
#refstarRA_deg=np.array([0.]*n_idx)
#refstarDEC_deg=np.array([0.]*n_idx)
#refstar_position=np.array([0.]*n_idx)
    
label_source=[['']*n_refstar]*n_idx
#print(label_source)
Rmag_targets=np.array([0.]*n_idx)
Rmag_err_targets=np.array([0.]*n_idx)
SNR_targets=np.array([0.]*n_idx)
#fwhm=np.array([0.]*n_idx)
#d_rmag_fitting=np.zeros((n_idx,n_refstar-1))
#radec_pix=np.array([[0.,0.]]*n_refstar)
bkg_err=np.array([0.]*n_idx)
bkg_err_ratio=np.array([0.]*n_idx)
R_square_rmag=np.array([0.]*n_idx)

err_all=np.array([0.]*n_idx)
#err_fitting_ratio=np.array([0.]*n_idx)
#mag_instrument_sum_err=np.array([0.]*n_idx)
err_count_per_img=np.array([0.]*n_idx)
err_mag_fitting_per_img=np.array([0.]*n_idx)
err_mag_instrument_per_img=np.array([0.]*n_idx)
err_mag_total_per_img=np.array([0.]*n_idx)


#sys.exit(0)

aper_sum0=np.zeros((n_idx,n_refstar))
aper_sum1=np.zeros((n_idx,n_refstar))
aper_sum_err0=np.zeros((n_idx,n_refstar))
aper_sum_err1=np.zeros((n_idx,n_refstar))
bg_subtracted_sum=np.zeros((n_idx,n_refstar))
bg_subtracted_err=np.zeros((n_idx,n_refstar))
bg_subtracted_sum_err=np.zeros((n_idx,n_refstar))
#bg_subtracted_sum_err_p=np.zeros((n_idx,n_refstar))
#bg_subtracted_sum_err_m=np.zeros((n_idx,n_refstar))
rmag_fitting=np.zeros((n_idx,n_refstar))
mag_instrument=np.zeros((n_idx,n_refstar))
mag_instrument_sum_err=np.zeros((n_idx,n_refstar))
#mag_instrument_sum_err_p=np.zeros((n_idx,n_refstar))
#mag_instrument_sum_err_m=np.zeros((n_idx,n_refstar))
mag_instrument_err=np.zeros((n_idx,n_refstar))
#mag_instrument_err_p=np.zeros((n_idx,n_refstar))
#mag_instrument_err_m=np.zeros((n_idx,n_refstar))

calendar_date=['']*n_idx


JD=df_info['JD'][idx_fitsheader]
ID=df_info['ID'][idx_fitsheader]
#JD=df_info['JD'][idx_fitsheader]
ID.tolist()

#sys.exit(0)
#=======================

j1=0
for j2 in idx_refstar:

#    print('j1,j2',j1,j2)
    refstarID[j1]=df_refstar['RefStarID'][j2]
    ra_deg[j1]=df_refstar['RefStarRA_deg'][j2]
    dec_deg[j1]=df_refstar['RefStarDEC_deg'][j2]
#    print(ra_deg[j1],dec_deg[j1])
    radec_deg[j1]=[ra_deg[j1],dec_deg[j1]]
#    print(j1,j2,radec_deg[j1])
    ra_hhmmss[j1]=df_refstar['RefStarRA_hhmmss'][j2]
    dec_ddmmss[j1]=df_refstar['RefStarDEC_ddmmss'][j2]
#    print(i,j1,radec_deg[j1],ra_hhmmss[j1],dec_ddmmss[j1])
    rmag[j1]=df_refstar['Rmag'][j2]
    rmag_err[j1]=df_refstar['Rmag_err'][j2]
    r_circle[j1]=df_refstar['R_circle_as'][j2]
    r_inner[j1]=df_refstar['R_inner_as'][j2]
    r_outer[j1]=df_refstar['R_outer_as'][j2]
#    position= SkyCoord(ra_deg[j1],dec_deg[j1],unit=(u.hourangle,u.deg),frame='icrs')
    position= SkyCoord(ra_deg[j1],dec_deg[j1],unit=(u.deg),frame='icrs')    
#    print(position)
    r_circle_as=r_circle[j1]*u.arcsec
#    print(r_circle_as)
    aperture[j1]=SkyCircularAperture(position, r=r_circle_as)
#    print(aperture[j1])
    r_inner_as=r_inner[j1]*u.arcsec
    r_outer_as=r_outer[j1]*u.arcsec
    aper_annu[j1]=SkyCircularAnnulus(position,r_inner_as,r_outer_as)
#    print(aper_annu[j1])
#    print(aper_annu[j1].r_in)
#    print(aper_annu[j1].r_out)

#    refstarID[j1]=df_refstar['RefStarID'][j2]
#    refstarRA_deg[j1]=df_refstar['RefStarRA_deg'][j2]
#    refstarDEC_deg[j1]=df_refstar['RefStarDEC_deg'][j2]
#    print(refstarID[j1])
#    print(refstarRA_deg[j1]
#    print(refstarDEC_deg[j1])
    
    j1=j1+1


dx=500
dy=dx

#=======================

k=0
for i in idx_fitsheader_canfit:
    print('=================================')   
    print('idx',i, ') ID =',ID[i],', #', k+1,'/',n_idx,'image of',obj_name)
    fits_root=fits_ori[i].split('.',-1)[0].split('_calib',-1)[0]
    fits_calib=fits_root+'_calib.fits'
    print('... filename :',fits_calib)
    #print(fits_root)
    #print(fits_calib)
    #print(fits_ori)

#   sys.exit(0)
#    print(radec_deg)
#    print(rmag)
    date=fits_root.split('@',-1)[0].split('-',-1)[-1]
    year=date[0:4]
    month=date[4:6]
    day=date[6:8]
    yearmonth=date[0:6]
#   sys.exit(0)
    dir_file=yearmonth+'/slt'+date+'_calib_sci/'
#    dir_reg=yearmonth+'/slt'+date+'_reg/'
    hdu=fits.open(dir_file+fits_calib)[0]
    imhead=hdu.header
    imdata=hdu.data    
    wcs = WCS(imhead)
#    print(wcs)

    calendar_date[k]=julian.from_jd(JD[i],fmt='jd')
    print('calendar_date',calendar_date[k])


#    radec_pix=wcs.all_world2pix(ra_deg,dec_deg,1)
    ra_pix,dec_pix=wcs.all_world2pix(ra_deg,dec_deg,1)
    ra_pix=ra_pix.tolist()
    dec_pix=dec_pix.tolist()
#    print(ra_pix,dec_pix)
#    print()


#    mask=np.array

#    background level and error
#    using FWHM
#    sigma_fwhm=fwhm[k]/2.35
#    print(k,fwhm[k],sigma_fwhm)


#    using the median absolute deviation to estimate the background noise level gives a value
#    bkg_noise=mad_std(imdata)
#    print(bkg_noise)
    print()
#    This method provides a better estimate of the background and background noise levels:
    bkg_mean,bkg_median,bkg_std=sigma_clipped_stats(imdata,sigma=3.)
    print('... bkg_mean,bkg_median,bkg_std')
    print('...','%.4f' %bkg_mean,'%.4f' %bkg_median,'%.4f' %bkg_std)
    # bkg_std = sqrt((bkg-bkg_mean)^2/n) = background noise
    
    imdata[imdata<bkg_median*0.9]=bkg_median
    
    bkg_err_ratio[k]=bkg_std/bkg_median
    print('... bkg_err_ratio = bkg_std/bkg_median =','%.4f' %bkg_err_ratio[k])    
#    bkg_err[k]=-2.5*np.log10(bkg_std/bkg_median)
#    print('bkg_err (mag)',bkg_err[k])  
#    1176.6646745343921 1176.0827349796416 39.1334108277639
#    np.argwhere(np.isnan(imdata))
#    time.sleep(1)
#    print(WCS.world_axis_physical_types)
    err_imdata=bkg_err_ratio[k]*imdata
    print()
    

    for j1 in range(n_refstar):
        print('idx',i, '-',j1+1,'/',n_refstar,', Ref.Star:',refstarID[j1],') #', k+1,'/',n_idx,'image of',obj_name )
#        aper_annu_pix[j1]=CircularAnnulus(ra_pix[j1],dec_pix[j1],r_inner_as,r_outer_as)
        aperture_pix[j1]=aperture[j1].to_pixel(wcs)
#        phot_table = aperture_photometry(imdata, aperture_pix[i],error=err_imdata)        
        aper_annu_pix[j1]=aper_annu[j1].to_pixel(wcs)
#        print(aper_annu_pix[j1])
        r_pix=aperture_pix[j1].r        
        r_in_pix=aper_annu_pix[j1].r_in
        r_out_pix=aper_annu_pix[j1].r_out
#        print(r_in_pix,r_out_pix)
#        position_pix=(ra_pix[j1],dec_pix[j1])
#        position_pix=np.transpose(positions_pix)
#        r_pix=aperture_pix[i].r        
#       print('r_pix =',r_pix)
        apper=[aperture_pix[j1],aper_annu_pix[j1]]
#        phot_annu_table = aperture_photometry(imdata, apper)
        phot_annu_table = aperture_photometry(imdata, apper,error=err_imdata)
#        phot_annu_table = aperture_photometry(imdata, aperture_pix[i],error=err_imdata)
        aper_annu_sum0=phot_annu_table['aperture_sum_0']
        print('... aper_annu_sum0 =','%.2f' %aper_annu_sum0)
        aper_annu_sum1=phot_annu_table['aperture_sum_1']
        print('... aper_annu_sum1 =','%.2f' %aper_annu_sum1)
        aper_sum0[k][j1]=aper_annu_sum0
        aper_sum1[k][j1]=aper_annu_sum1

        aper_annu_sum_err0=phot_annu_table['aperture_sum_err_0']
        print('... aper_annu_sum_err0 =','%.2f' %aper_annu_sum_err0)
        aper_annu_sum_err1=phot_annu_table['aperture_sum_err_1']        
        print('... aper_annu_sum_err1 =','%.2f' %aper_annu_sum_err1)
        aper_sum_err0[k][j1]=aper_annu_sum_err0
        aper_sum_err1[k][j1]=aper_annu_sum_err1
        
        area_sum0=aperture_pix[j1].area
        print('... area of sum0 =','%.2f' %area_sum0)
        area_sum1=aper_annu_pix[j1].area
        print('... area of sum1 =','%.2f' %area_sum1)
        
#        area_annu=aper_annu_pix[j1].area
#        print('area_annu',area_annu)
#        bkg_mean = phot_annu_table['aperture_sum_1'] / area_annu
#        print('bkg_mean',bkg_mean)
#        bkg_sum = bkg_mean * aperture_pix[j1].area
        
        annulus_mask = aper_annu_pix[j1].to_mask(method='center')
#        plt.imshow(annulus_mask)
#        plt.colorbar()
#        bkg_median = []
#        for mask in annulus_masks:
        annulus_data = annulus_mask.multiply(imdata)
        mask = annulus_mask.data
        annulus_data_1d = annulus_data[mask > 0]
        _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
        print('... median_sigclip (annual mask per pix)  =','%.2f' %median_sigclip)
        print('... bkg_mean (imdata statistic per pix)   =','%.2f' %bkg_mean)
        print('... bkg_median (imdata statistic per pix) =','%.2f' %bkg_median)
        bkg_sum=median_sigclip*area_sum0
        print('... bkg_sum (median_sigclip*area_sum0)           =','%.2f' %bkg_sum)
#        bkg_median.append(median_sigclip)
#        bkg_median = np.array(bkg_median)
        bkg_sum=aper_annu_sum1/area_sum1*area_sum0
        print('... bkg_sum (aper_annu_sum1/area_sum1*area_sum0) =','%.2f' %bkg_sum)
#        print()
        
#        phot_annu_table['bkg_sum']=bkg_sum
        bg_subtracted_sum[k][j1]=aper_annu_sum0-bkg_sum        
        print('... bg_subtracted_sum =','%.2f' %bg_subtracted_sum[k][j1])

        bg_subtracted_err[k][j1]=np.sqrt(aper_annu_sum_err0**2+aper_annu_sum_err1**2)
        print('... bg_subtracted_err =','%.2f' %bg_subtracted_err[k][j1])
        
#        print(phot_annu_table.colnames)
        phot_annu_table['bg_subtracted_sum'] = bg_subtracted_sum[k][j1]
        phot_annu_table['bg_subtracted_err'] = bg_subtracted_err[k][j1]
#        print()
        
#        phot_annu_table['bg_subtracted_sum[k][j1]_err'] = bg_subtracted_sum[k][j1]
#        print(phot_annu_table)        
        mag_instrument[k][j1]=-2.5*np.log10(bg_subtracted_sum[k][j1])
#        if bg_subtracted_sum[k][j1]>0:            
#            mag_instrument[k][j1]=-2.5*np.log10(bg_subtracted_sum[k][j1])
#        else:
#            rmag_uplims[j1]=True
#            mag_instrument[k][j1]=-2.5*np.log10(bkg_sum)
        print('... mag_instrument (sum) =','%.4f' %mag_instrument[k][j1])        
        bg_subtracted_sum_err[k][j1]=bg_subtracted_sum[k][j1]+bg_subtracted_err[k][j1]
        mag_instrument_sum_err[k][j1]=-2.5*np.log10(bg_subtracted_sum_err[k][j1])
        print('... mag_instrument (sum+err) =','%.4f' %mag_instrument_sum_err[k][j1])
        mag_instrument_err[k][j1]=abs(mag_instrument_sum_err[k][j1]-mag_instrument[k][j1])
        print('... mag_instrument (err) =','%.4f' %mag_instrument_err[k][j1])

#        bg_subtracted_sum_err_p[k][j1]=bg_subtracted_sum[k][j1]+bg_subtracted_err[k][j1]
#        mag_instrument_sum_err_p[k][j1]=-2.5*np.log10(bg_subtracted_sum_err_p[k][j1])
#        print('... mag_instrument_sum+err =','%.4f' %mag_instrument_sum_err_p[k][j1])
#        mag_instrument_err_p[k][j1]=abs(mag_instrument_sum_err_p[k][j1]-mag_instrument[k][j1])
#        print('... mag_instrument+err =','%.4f' %mag_instrument_err_p[k][j1])
        
#        bg_subtracted_sum_err_m[k][j1]=bg_subtracted_sum[k][j1]-bg_subtracted_err[k][j1]
#        mag_instrument_sum_err_m[k][j1]=-2.5*np.log10(bg_subtracted_sum_err_m[k][j1])
#        print('... mag_instrument_sum-err =','%.4f' %mag_instrument_sum_err_m[k][j1])
#        mag_instrument_err_m[k][j1]=abs(mag_instrument_sum_err_m[k][j1]-mag_instrument[k][j1])
#        print('... mag_instrument-err =','%.4f' %mag_instrument_err_p[k][j1])        

#        aper_sum_err0[k][j1]=phot_annu_table['aperture_sum_err_0']
#        aper_sum_err1[k][j1]=phot_annu_table['aperture_sum_err_1']
        print()
            
        for col in phot_annu_table.colnames:
            phot_annu_table[col].info.format = '%.8g'  # for consistent table output
#            phot_annu_table[col].info.format = '%.2f'  # for consistent table output
        print(phot_annu_table)
        print()
#    print(phot_annu_table.colnames)
    print('=================================')

#    bg_subtracted_sum_mean,bg_subtracted_sum_median,bg_subtracted_sum_std=sigma_clipped_stats(bg_subtracted_sum,sigma=3.)
#    print('bg_subtracted_sum_mean,bg_subtracted_sum_median,bg_subtracted_sum_std')
#    print(bg_subtracted_sum_mean,bg_subtracted_sum_median,bg_subtracted_sum_std)
    

#    print('bg_subtracted_err[k]')
#    print(bg_subtracted_err[k])
    err_count_per_img[k]=np.sqrt(sum(bg_subtracted_err[k]**2)/n_refstar)
    print('... measurement error after background subtraction (counts):','%.4f' %err_count_per_img[k])  
    err_mag_instrument_per_img[k]=np.sqrt(sum(mag_instrument_err[k]**2)/n_refstar)
    print('... magnitude error after background subtraction (mag):','%.4f' %err_mag_instrument_per_img[k])  


#    print('aperture sum',aper_sum0[k])
#    print('aperture sum error',aper_sum_err0[k])
#    print('annulus sum',aper_sum1[k])    
#    print('annulus sum error',aper_sum_err1[k])


#    mag_instrument[k][1:]=mag_instrument[k][1:]
#    print(mag_instrument[k][1:])
#    print(rmag[1:])

#    mag_instrument_source=mag_instrument
#    rmag_source=rmag

#    no values    
#    mag_instrument_bkg=mag_instrument[-1]
#    rmag_bkg=rmag[-1]
    
    
    

    '''  
#    aperture_pix[i]=aperture.to_pixel(wcs)
#    print(aperture_pix[i])
#    positions_pix=aperture_pix[i].positions
#    r_pix=aperture_pix[i].r
#    print('r_pix =',r_pix)
#    positions_pix=(ra_pix,dec_pix)
#    positions_pix=np.transpose(positions_pix)
#    aperture_pix[i]=CircularAperture(positions_pix, r_pix)
    
#    mask_circle=aperture_pix[i].to_mask(method='center')
#    mask_data=mask_circle[0].multiply(imdata)
    
    
    err_imdata=bkg_err_ratio[k]*imdata
#    phot_table = aperture_photometry(imdata, aperture,wcs=wcs)
    phot_table = aperture_photometry(imdata, aperture_pix[i],error=err_imdata)

    '''

    '''
#    aper_annu_pix[i]=aper_annu.to_pixel(wcs)
#    print(aper_annu_pix[i])
#    r_in_annu_pix=aper_annu_pix[i].r_in
#    r_out_annu_pix=aper_annu_pix[i].r_out
#    print(r_in_annu_pix,r_out_annu_pix)

   '''  
 
    '''
    annulus_masks = aper_annu_pix[i].to_mask(method='center')
    bkg_median = []
    for mask in annulus_masks:
        annulus_data = mask.multiply(imdata)
        annulus_data_1d = annulus_data[mask.data > 0]
        _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
        bkg_median.append(median_sigclip)
    bkg_median = np.array(bkg_median)
    '''
#    bkg_sum = bkg_median * aperture_pix[i].area
#    bkg_sum = bkg_mean * aperture_pix[i].area
#    print(bkg_sum)
#    phot_table['bkg_sum']=bkg_sum
#    phot_table['bkg_sum'].info.format = '%.8g'  # for consistent table output

#    final_sum = phot_annu_table['aperture_sum_0'] - bkg_sum
#    final_sum=aper_sum-bkg_sum
#   print(final_sum)
#    phot_annu_table['residual_aperture_sum'] = final_sum
#    phot_table['bg_subtracted_sum'] = final_sum  
#    phot_annu_table['residual_aperture_sum'].info.format = '%.8g'  # for consistent table output
#   print(phot_annu_table['residual_aperture_sum'])  
#    print(phot_annu_table)

#    phot_table['bg_subtracted_sum'].info.format = '%.8g'  # for consistent table output
#    print(phot_table)

#   sys.exit(0)
    '''
    j=0
    rmag[0]=-9999.
    rmag_fitting[0]=rmag[0]
    print()
    print('j final_sum         mag_instrument      Rmag')
    for j in range(n_refstar):
        mag_instrument[j]=-2.5*np.log10(bg_subtract_sum[j])
        print(j, final_sum[j],mag_instrument[j],rmag[j])

    mag_instrument[k][1:]=mag_instrument[1:]
    rmag[1:]=rmag[1:]

    mag_instrument_source=mag_instrument
    rmag_source=rmag
    '''
#    no values    
#    mag_instrument_bkg=mag_instrument[-1]
#    rmag_bkg=rmag[-1]

#    print('Instrument Mag of RefStars',mag_instrument[k][1:])
#    print('Rmag of RefStars',rmag[1:])    
    print()
#    shift,slope=polyfit(mag_instrument[k][1:],rmag[1:],1)

#    param,res=polyfit(mag_instrument[k][1:],rmag[1:],1,full=True)
    shift,slope=polyfit(mag_instrument[k][1:],rmag[1:],1,w=1/(mag_instrument_err[k][1:]))
#    shift,slope=polyfit(mag_instrument[k][1:],rmag[1:],1)

#    residual=res[0].tolist()[0]
#    print('shift =','%.3f' %shift)
#    print('slope =','%.3f' %slope)
#    print('residual =', residual)

    rmag_fitting[k]=shift+slope*mag_instrument[k]
    print('... rmag_fitting (Ref.Stars only):',rmag_fitting[k][1:])
    print('... rmag_fitting (Target only):',rmag_fitting[k][0])
    
    Rmag_targets[k]=rmag_fitting[k][0]

#    SN=signal_mean/signal_std
#    SN=bkg_mean/bkg_std ??    
#    ErrorRmag=2.5*log10(1+N/S)=1.086/(S/N)
#    signal=final_sum
#    noise=bkg_std*aperture_pix[i].area
#    mag_err=1.086/snr

    
#    Rmag_err_targets[k]=mag_err[0]

    '''
#    plt.figure(figsize=(8,6))
    mask_circle=aperture_pix[j1].to_mask(method='center')
#    plt.imshow(mask_circle[3])
#    plt.colorbar()
#    plt.show()
    data_in_mask=mask_circle[3].multiply(imdata)
#    plt.imshow(data_in_mask)
#    plt.colorbar()
#    plt.show()
#    plt.close()
    
    
    mask_median = []
    for mask in mask_circle:
        circle_data = mask.multiply(imdata)
        circle_data_1d = circle_data[mask.data > 0]
        _, median_sigclip, _ = sigma_clipped_stats(circle_data_1d)
        mask_median.append(median_sigclip)
#        bkg_median = np.array(bkg_median)
    '''

    

#    snr=mask_median/bkg_std
#    snr=bg_subtracted_sum/aperture_pix[j1].area/bkg_std
#    snr=10
#    mag_err=1.086/snr
    
#    print('snr',snr)
#    print('mag_err',mag_err)
#    mag_err_target=mag_err
#    Rmag_err_targets[k]=mag_err
#    SNR_targets[k]=snr 

    '''
    phot = aperture_photometry(imdata, aperture_pix[i])
#    phot['circle_median'] = bkg_median
#    phot['bkg_median'] = bkg_median
    phot['aper_bkg'] = bkg_median * aperture_pix[i].area
#    phot['aper_bkg'] = bkg_sum
    phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']
    for col in phot.colnames:
        phot[col].info.format = '%.8g'  # for consistent table output
    print(phot)
    ''' 


#    refstar_ra_pix,refstar_dec_pix=wcs.all_world2pix(refstarRA_deg,refstarDEC_deg,1)
#    refstar_ra_pix=refstar_ra_pix.tolist()
#    refstar_dec_pix=refstar_dec_pix.tolist()
 
    
    
    figtitle='(ID='+str(ID[i])+', JD='+str('%.0f' %JD[i])+') '+fits_ori[i]

            
    fig1=plt.figure(figsize=(8,4))
    fig1.suptitle(figtitle,fontsize=font_middle)
#    fig1 = plt.figure(constrained_layout=True)
#    widths = [1,1]
#    heights = [1]
#    spec1 = fig1.add_gridspec(ncols=len(widths), nrows=len(heights), width_ratios=widths,
#                          height_ratios=heights)
    
#    fig1,axs=plt.subplots(1,2)
#    axs=fig1.add_subplot(1,2,1)
    plt.subplot(121)
#    plt.tight_layout(rect=[0,0,1,1])

#    plt.figure(dpi=200)
#    if change ra,dec of reference star, run 'gasp_refStar_radec.py again.    
    norm = simple_norm(imdata, 'sqrt', percent=99)
    plt.imshow(imdata, norm=norm)
#    im=axs[0].imshow(imdata,original='lower',norm=norm)
#    fig1.colorbar(im)
#    axs[0].plot(imdata)
    
    for i in range(n_refstar):
        aper_annu_pix[i].plot(color='black', lw=1)      
        if i==0:
            aperture_pix[i].plot(color='red', lw=1)
            plt.text(ra_pix[i]*1.02,dec_pix[i]*0.98,str(refstarID[i]),fontsize=font_middle,color='red')
        else:
            aperture_pix[i].plot(color='white', lw=1)
            plt.text(ra_pix[i]*1.02,dec_pix[i]*0.98,str(refstarID[i]),fontsize=font_middle)
    

    x0=ra_pix[0]
    y0=dec_pix[0]
    x1=x0-dx
    x2=x0+dx
    y1=y0-dy
    y2=y0+dy
    plt.xlim(x1,x2)
    plt.ylim(y1,y2)
    plt.xlabel('RA (J2000)',fontsize=font_small)
    plt.ylabel('Dec (J2000)',fontsize=font_small)
    plt.xticks(fontsize=font_small)
    plt.yticks(fontsize=font_small)    
#    plt.title(figtitle,fontsize=font_small)
    text_rmag='Rmag='+str('%.3f' %rmag_fitting[k][0])
#    text_rmag_err='+/-'+str('%.3f' %mag_err_target)
#    text_fwhm='fwhm='+str('%.3f' %fwhm[i])
#    text_radius='r = fwhm x '+str(xfwhm)
#    text_fig=text_rmag+', '+text_fwhm
#    plt.text(x0*0.95,y0*0.94,obj_name,fontsize=font_small,color='red')
#    plt.text(x0*0.95,y0*0.90,text_rmag,fontsize=font_small)
#    for i in range(n_refstar):
#        plt.text(ra_pix[i]*0.95,dec_pix[i]*0.95,str(refstarID[i]),fontsize=font_small)

#    plt.text(x0*0.95,y0*0.90,'refstar',fontsize=font_small)
#    plt.text(x0*0.95,y0*0.93,text_rmag+text_rmag_err,fontsize=font_small)
#    plt.text(x0*0.95,y0*0.90,text_fwhm,fontsize=font_small)
#    plt.text(x0*0.95,y0*0.87,text_radius,fontsize=font_small)
#    plt.text(x2*0.55,y2*0.45,text_fwhm,fontsize=font_small)
#    plt.show()    
#    plt.savefig(dir_obj+'Rmag_InstMag_aper_'+obj_name+'_'+date+'_fwhm'+str(xfwhm)+'_'+str(k)+'.png')
#    plt.savefig(dir_obj+'Rmag_InstMag_aper_'+obj_name+'_annu_'+date+'_'+str(k)+'.png')
#    plt.close()


    plt.subplot(122,aspect=0.5)
#    plt.subplots_adjust(top=0.75)
#    plt.tight_layout(rect=[0,0,0.75,0.5])
#    plt.figure()
#    plt.figure(figsize=(8,6))
#    plt.scatter(mag_instrument[k][1:],rmag[1:])
#    plt.plot(mag_instrument[k][1:],rmag[1:],'o')

    d_rmag_fitting=rmag[1:]-rmag_fitting[k][1:]
    print('... fitting residual:',d_rmag_fitting)

#    rmag=shift+slope*mag_instrument
#    mag_instrument_fitting_refstar=(rmag_fitting[k][1:]-shift)/slope
#    count_instrument_fitting_refstar=10**(mag_instrument_fitting_refstar/(-2.5))
#    d_count_instrument_fitting=abs(count_instrument_fitting_refstar-bg_subtracted_sum[1:])
#    d_mag_instrument_fitting=-2.5*np.log10(d_count_instrument_fitting)
    err_mag_fitting_per_img[k]=np.sqrt(sum(d_rmag_fitting**2)/n_refstar)    
    print('... fitting error (mag) =','%.4f' %err_mag_fitting_per_img[k])
#    print('measurement error (mag)',bkg_err[k])    
#    SSregression_rmag=sum((rmag_fitting[k][1:]-np.mean(rmag[1:]))**2)
    SSE_rmag=sum((d_rmag_fitting)**2)
    SST_rmag=sum((rmag[1:]-np.mean(rmag[1:]))**2)
    R_square_rmag[k]=1-SSE_rmag/SST_rmag
    print('... error of measurement (ratio of bkg_std/bkg_mediean) =','%.4f' %bkg_err_ratio[k])
    err_mag_total_per_img[k]=np.sqrt((err_mag_fitting_per_img[k])**2+(err_mag_instrument_per_img[k])**2)
    print('... error of (calibration + observation) [mag]:','%.4f' %err_mag_total_per_img[k])
    text_fitting='Rmag='+str('%.2f' %shift)+'+'+str('%.2f' %slope)+'*(Inst.Mag), $R^2$='+str('%.4f' %R_square_rmag[k])
    label_target=obj_name+', Rmag='+str('%.2f' %rmag_fitting[k][0])+'$\pm$'+str('%.2f' %err_mag_total_per_img[k])
    for i in range(n_refstar):
        label_source[k][i]=obj_name+', Rmag='+str('%.2f' %rmag_fitting[k][i])
    print('label_source[k]',label_source[k])
#    label_source[k]=obj_name+', Rmag='+str('%.2f' %rmag_fitting[k])
    plt.plot(mag_instrument[k],rmag_fitting[k],'-',color='grey',label=text_fitting)
    plt.errorbar(mag_instrument[k][1:],rmag[1:],yerr=mag_instrument_err[k][1:],fmt='o',label='Reference Stars')     
#    plt.plot(mag_instrument[k][0],rmag_fitting[k][0],'o',color='red',label=label_target)
#    plt.errorbar(mag_instrument[k][0],rmag_fitting[k][0],yerr=err_mag_fitting_per_img[k],fmt='o',color='red',label=label_target)
    plt.errorbar(mag_instrument[k][0],rmag_fitting[k][0],yerr=mag_instrument_err[k][0],fmt='o',color='red',label=label_target)

    print('... bg_subtracted_sum =',bg_subtracted_sum[k][j1])
#    if bg_subtracted_sum[k][j1]>0:
#        plt.errorbar(mag_instrument[k][0],rmag_fitting[k][0],yerr=mag_instrument_err[k][0],fmt='o',color='red',label=label_target)
#    else:
#        plt.errorbar(mag_instrument[k][0],rmag_fitting[k][0],yerr=0.2,uplims=True,fmt='o',color='red',label=label_target)        
#       plt.plot(mag_instrument[k][1:],rmag[1:],'o',label='Ref.Star')
#        plt.plot(mag_instrument[k][0],rmag_fitting[k][0],'o',color='red',label=label_target)
#    else:
#        plt.errorbar(mag_instrument[k][0],rmag_fitting[k][0],yerr=0.2,uplims=['True'],marker='o',color='red',label=label_target)
#    plt.plot(mag_instrument[k][1:],rmag_fitting[k][1:],'o',color='steelblue',label=label_source[1:])
#    plt.errorbar(mag_instrument[k][0],rmag_fitting[k][0],yerr=mag_instrument_err[k][0],fmt='o',color='red',label=label_target)
#    plt.errorbar(mag_instrument_source,rmag_source,yerr=mag_err[:-1])


    xx=mag_instrument[k]
    yy=rmag_fitting[k]
    xxmin=min(xx)
    xxmax=max(xx)
    yymin=min(yy)
    yymax=max(yy)
    xxwidth=abs(xxmin-xxmax)
    yywidth=abs(yymin-yymax)
#    xxrange=range(int(xxmin*10),int(xxmax*10),5)
    print('... InstMag_width =','%.2f' %xxwidth,', Rmag_width =','%.2f' %yywidth)
    
    for i in range(n_refstar):
        if i==0:
            plt.text((xx[i]-0.01*xxwidth),(yy[i]+0.05*yywidth),str(refstarID[i]),fontsize=font_small,color='red')
        else:
            plt.text((xx[i]-0.01*xxwidth),(rmag[i]),str(refstarID[i]),fontsize=font_small,color='steelblue')
    plt.xlabel('Instrument Mag',fontsize=font_small)
    plt.ylabel('Rmag',fontsize=font_small)
#    plt.title(figtitle,fontsize=font_small)
#    plt.xticks(np.arange(min(xx)-1,max(xx)+1,0.5),fontsize=font_small)
#    plt.yticks(np.arange(min(yy)-1,max(yy)+1,0.5),fontsize=font_small)
    plt.xticks(fontsize=font_small)
    plt.yticks(fontsize=font_small)    
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()


    

#    print(mag_instrument,Rmag_targets)
#    plt.autoscale(enable=True,axis='both')
#    ax=plt.gca()
#    ax.relim()
#    ax.autoscale()
#    ax.Axes.autoscale #(enable=True,axis='both')
    
#    set_xlim(auto=True)
#    set_ylim(auto=True)
#    plt.show()
#    plt.legend(loc='best') #,fontsize=font_small)

    plt.legend(loc='lower left',bbox_to_anchor=(0.0,1.02),fontsize=font_small)
#    plt.legend(loc='upper left', bbox_to_anchor=(-0.,1.65),fontsize=font_small)
#    plt.savefig(dir_obj+'Rmag_InsMag_rmag_'+obj_name+'_annu_'+date+'_'+str(k)+'.png')
    plt.savefig(dir_obj+'Rmag_InsMag_'+obj_name+'_annu_'+date+'_'+str(k+1)+'.png',dpi=300)    
    plt.close()    
#    print('target new Rmag =',rmag_fitting[k][0])
#    print('Instrument Mag',mag_instrument)
#    print('=================================')
    k=k+1
    


#print('=================================')
print()    
print('Rmag',Rmag_targets)
print()   

df_out=df_info.iloc[idx_fitsheader_canfit]
pd.options.mode.chained_assignment = None  # default='warn'
#df_out['Rmag']=pd.Series(Rmag_targets,index=df_out.index)
df_out['Rmag']=Rmag_targets
df_out['Rmag']=df_out['Rmag'].map('{:.4f}'.format)
df_out['ErrorCounts']=err_count_per_img
df_out['ErrorCounts']=df_out['ErrorCounts'].map('{:.4f}'.format)
df_out['ErrorInstMag']=err_mag_instrument_per_img
df_out['ErrorInstMag']=df_out['ErrorInstMag'].map('{:.4f}'.format)
df_out['ErrorFitting']=err_mag_fitting_per_img
df_out['ErrorFitting']=df_out['ErrorFitting'].map('{:.4f}'.format)
df_out['ErrorRmag']=err_mag_total_per_img
df_out['ErrorRmag']=df_out['ErrorRmag'].map('{:.4f}'.format)


#fmt="{:,.4f}"  # 1,234
#fmt="{:.4f}"  # 1234
#df_out.style.format({'RA_deg':fmt,'DEC_deg':fmt,'RA_pix':fmt,'DEC_pix':fmt,'Zmag':fmt,'FWHM':fmt,'Altitude':fmt,'Airmass':fmt,'Rmag':fmt,'ErrorCounts':fmt,'ErrorInstMag':fmt,'ErrorFitting':fmt,'ErrorMagTotal':fmt})
#df_out.style.format({'Rmag':fmt,'ErrorCounts':fmt,'ErrorInstMag':fmt,'ErrorFitting':fmt,'ErrorMagTotal':fmt})

print(df_out)
file_Rmag_out_full=dir_obj+'Rmag_aperture_'+obj_name+'_annu.txt'
df_out.to_csv(file_Rmag_out_full,sep='|')
# -----------------------------------
df_Rmag=df_out[['JD','Rmag','ErrorRmag']].reset_index(drop=True)
df_Rmag['JD']=df_Rmag['JD'].map('{:.4f}'.format)
file_Rmag_out_simple=dir_obj+'Rmag_'+obj_name+'_all.txt'
df_Rmag.to_csv(file_Rmag_out_simple,sep=',',float_format='%.4f')
# -----------------------------------
#=======================

#df_baddata_note=pd.read_csv('bad_data_note.txt',sep='\t')
#print(df_baddata_note)

df_baddata=df_baddata_note.loc[(df_baddata_note['NoteIdx']==2)].reset_index(drop=True)
n_baddata=len(df_baddata)
print('... there are',n_baddata,'bad fits')
list_baddata=df_baddata['filename'].tolist()
#print(list_baddata)
#sys.exit(0)     

n_idx=len(idx_fitsheader_canfit)

idx_skip=[]
print('before drop bad fits: idx_fitsheader_canfit =',idx_fitsheader_canfit)
for i in range(n_idx):
    j=idx_fitsheader_canfit[i]
    fitsname=fits_ori[j]
#    print(j,fitsname)
    if fitsname in list_baddata:
        idx_skip.append(i)
#        print(idx_skip)
        #ignore=ignore+1
        print(i,j,fitsname,'... bad fits')
#        print('skip',j,fitsname,'... bad fits ...')
        print()
    else:
        print(i,j,fitsname)


#del idx_fitsheader_canfit[idx_skip]
idx_fitsheader_keep=np.delete(idx_fitsheader_canfit,idx_skip)
print(idx_fitsheader_keep)
ID_keep=idx_fitsheader_keep+1
print(ID_keep)
idx_fitsheader_skip=idx_fitsheader_canfit[idx_skip]
ID_skip=idx_fitsheader_skip+1
print(ID_skip)
n_idx_keep=len(idx_fitsheader_keep)
n_skip=len(idx_skip)
#JD_keep=JD
print()
print('... skip idx_fitsheader_keep',idx_skip,'=',idx_fitsheader_canfit[idx_skip],'... bad fits ...')
print()


print('after drop bad fits: idx_fitsheader_keep =',idx_fitsheader_keep)
print('... keep fits ...')
print(fits_ori[idx_fitsheader_keep])    
        
#=======================
#sys.exit(0)
df_out_keep=df_out
df_out_keep=df_out_keep[~df_out_keep['ID'].isin(ID_skip)]
print(df_out_keep)
# -----------------------------------
df_Rmag_keep=df_out_keep[['JD','Rmag','ErrorRmag']].astype(float).reset_index(drop=True)
#df_Rmag_keep['JD']=df_Rmag_keep['JD'].map('{:.4f}'.format)

print(df_Rmag_keep)
file_Rmag_out_keep=dir_obj+'Rmag_'+obj_name+'.txt'
df_Rmag_keep.to_csv(file_Rmag_out_keep,sep=',',index=False,float_format='%.4f')
# -----------------------------------
#JD_keep=df_out_keep['JD']
#Rmag_keep=df_out_keep['Rmag']
#ErrorRmag=df_out_keep['ErrorRmag']
# -----------------------------------
#sys.exit(0)


# -----------------------------------
fig2=plt.figure(figsize=(6,8))
#fig,axs=plt.subplots(2,1,sharex=True) #,gridspec_kw={'hspace':0})
fig2.suptitle(obj_name)
#plt.subplot(122,aspect=0.5)
#plt.subplots(nrows=3,ncols=1,sharex=True)
#ax1=plt.subplot(3,1,1)
#ax2=plt.subplot(311)
#f,(ax1,ax2,ax3)=plt.subplots(3,1,sharex=True,figsize=(6,8),constrained_layout=True)
plt.subplot(311)
plt.subplots_adjust(right=0.75)


#plt.figure()
#plt.scatter(JD,Rmag_targets,label='data')


for i in range(n_refstar):
    print(i)
    print(JD.tolist())
    print(rmag_fitting[:,i]) 
#    plt.scatter(JD,rmag_fitting[:,i]) # ,fmt='o')
#    plt.scatter(JD,rmag_fitting[:,i],label=str(refstarID[i]),marker='o')
#    plt.plot(JD,rmag_fitting[:,i],label=str(refstarID[i]),fmt='o')
#    plt.plot(JD,rmag_fitting[:,i],label=label_source[:,i],fmt='o')
    if i==0:
        print('i=',i)
        plt.errorbar(JD,rmag_fitting[:,i],yerr=err_mag_total_per_img,label=str(refstarID[i]),marker='.',color='red',linestyle='none')
#        plt.errorbar(JD,rmag_fitting[:,i],yerr=err_mag_total_per_img,label=str(refstarID[i]),marker='o',color='red',linestyle='none',fmt='+', mfc='white')
#	plt.errorbar(JD,Rmag_targets,yerr=err_mag_total_per_img,label=label_target,fmt='o')
    else:
        print('i=',i)
        plt.axhline(y=rmag[i],color=colors[i],lw=1)
#        plt.axhline(y=rmag[i]-0.1,color=colors[i],linestyle='--',lw=1)
#        plt.axhline(y=rmag[i]+0.1,color=colors[i],linestyle='--',lw=1)
        plt.scatter(JD,rmag_fitting[:,i],label=str(refstarID[i]),marker=markers[i],color=colors[i],facecolors='none')
        
        

#	print(i)
#	plt.scatter(JD,rmag_fitting[:,i],label=label_source[:,i],fmt='o')
#        plt.scatter(JD,rmag_fitting[:,i],label=str(refstarID[i]),marker='o')
#        plt.axhline(y=rmag[1:])
#	plt.axhline(JD,rmag[1:],linestyle='-')
#	for j in range(k):
#	    plt.scatter(JD,rmag_fitting[k][i],label=label_source[k][i],fmt='o')
#plt.plot(JD,Rmag_targets,'-')
        
#axs[0].set_xlabel('JD',fontsize=font_small)
plt.ylabel('Rmag',fontsize=font_small)
plt.gca().invert_yaxis()
#axs[0].set_title(obj_name) #,fontsize=font_small)
plt.xticks(fontsize=font_small)
plt.yticks(fontsize=font_small)
#tick.label.set_fontsize(14)
#axs[0].legend(loc='best') #,fontsize=font_small)

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#axs[0].legend(loc=7) 
#axs[0].tight_layout()
#axs[0].subplots_adjust(right=0.75) 
#plt.show()
#plt.savefig(dir_obj+'Rmag_JD_'+obj_name+'_RefStar_annu.png')
#plt.close()
#sys.exit(0)
# -----------------------------------
plt.subplot(312) #,sharex=ax1)
plt.errorbar(JD,Rmag_targets,yerr=err_mag_total_per_img,label=obj_name,fmt='.',color='red')

#plt.xlabel('JD',fontsize=font_small)
plt.ylabel('Rmag',fontsize=font_small)
plt.gca().invert_yaxis()
#axs[1].set_title(obj_name) #,fontsize=font_small)
plt.xticks(fontsize=font_small)
plt.yticks(fontsize=font_small)
#axs[1].legend(loc='best') #,fontsize=font_small)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# -----------------------------------

plt.subplot(313)
#plt.figure(figsize=(8,6))
#plt.figure()
#plt.plot(JD,Rmag_targets,'-',label='lightcurve',color='grey')
#plt.scatter(JD,Rmag_targets,label='data')
JD_keep=df_Rmag_keep['JD'] #.tolist()
Rmag_targets_keep=df_Rmag_keep['Rmag'] #.tolist()
err_mag_total_per_img_keep=df_Rmag_keep['ErrorRmag'] #.tolist()
#plt.scatter(JD_keep,Rmag_targets_keep,marker='.',color='red') #,label=obj_name,fmt='.',color='red')

plt.errorbar(JD_keep,Rmag_targets_keep,yerr=err_mag_total_per_img_keep,label=obj_name,fmt='.',color='red')


plt.xlabel('JD',fontsize=font_small)
plt.ylabel('Rmag',fontsize=font_small)
plt.gca().invert_yaxis()
#axs[1].set_title(obj_name) #,fontsize=font_small)
plt.xticks(fontsize=font_small)
plt.yticks(fontsize=font_small)
#axs[1].legend(loc='best') #,fontsize=font_small)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

#plt.show()
plt.savefig(dir_obj+'Rmag_JD_'+obj_name+'_RefStar_annu.png',dpi=300)
plt.close()
#sys.exit(0)

# -----------------------------------
'''
fig3=plt.figure
JD_keep=df_Rmag_keep['JD'] #.tolist()
Rmag_targets_keep=df_Rmag_keep['Rmag'] #.tolist()
err_mag_total_per_img_keep=df_Rmag_keep['ErrorRmag'] #.tolist()
#plt.scatter(JD_keep,Rmag_targets_keep,marker='.',color='red') #,label=obj_name,fmt='.',color='red')

plt.errorbar(JD_keep,Rmag_targets_keep,yerr=err_mag_total_per_img_keep,label=obj_name,fmt='.',color='red')


plt.xlabel('JD',fontsize=font_middle)
plt.ylabel('Rmag',fontsize=font_middle)
plt.gca().invert_yaxis()
#axs[1].set_title(obj_name) #,fontsize=font_small)
plt.xticks(fontsize=font_middle)
plt.yticks(fontsize=font_middle)
#axs[1].legend(loc='best') #,fontsize=font_small)
plt.legend(loc='best')

#plt.show()
plt.savefig(dir_obj+'Rmag_JD_'+obj_name+'.png',dpi=120)
plt.close()
'''
# -----------------------------------


rmag_nan=([])

k=0
for i in idx_fitsheader_keep:
    ID_fitsheader=ID[i]
    if Rmag_targets[k]==np.nan:
        rmag_nan=np.append(rmag_nan,ID_fitsheader)
    k=k+1

print()
print('Rmag with nan value:',rmag_nan)
print()
print('... write file to: ./'+file_Rmag_out_full)
print('... write file to: ./'+file_Rmag_out_simple)
print('... finished ...')
#==============================



'''
annulus_masks = aper_annu_pix[i].to_mask(method='center')
plt.imshow(annulus_masks)
plt.colorbar()

annulus_data = annulus_masks[0].multiply(imdata)
plt.imshow(annulus_data)
plt.colorbar()

mask = annulus_masks[0].data
annulus_data_1d = annulus_data[mask > 0]
annulus_data_1d.shape
'''



'''
annulus_aperture = CircularAnnulus(position, r_in=fwhm*6.* u.arcsec, r_out=fwhm*8.* u.arcsec)
apers=[aperture,annulus_aperture]
phot_table2=aperture_photometry(imdata,apers)
for col in phot_table2.colnames:
    phot_table2[col].info.format = '%.8g'  # for consistent table output
(phot_table2)

#bkg_mean = phot_table2['aperture_sum_1'] / annulus_aperture.area
'''
