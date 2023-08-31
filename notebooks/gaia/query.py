#!/usr/bin/env python
# encoding: utf-8

""" GAIA 2 mask """

import os
import sys
from contextlib import contextmanager

import numpy as np
import sep

from astropy import units as u
from astropy.table import Column
from astropy.units import Quantity
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

__all__ = ['image_gaia_stars', 'gaia_star_mask', 'gaia2mask']


@contextmanager
def suppress_stdout():
    """Suppress the output.

    Based on: https://thesmithfam.org/blog/2012/10/25/temporarily-suppress-console-output-in-python/
    """
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout
            
def image_gaia_stars(image, wcs, CustomQuery, radius=None, center=None, pixel=0.168,
                     mask_a=694.7, mask_b=4.04,
                     verbose=False, size_buffer=1.4,
                     tap_url=None):
    # Central coordinate
    if center is None:
        ra_cen, dec_cen = wcs.wcs_pix2world(
            image.shape[0] / 2, image.shape[1] / 2, 0)
        img_cen_ra_dec = SkyCoord(
            ra_cen, dec_cen, unit=('deg', 'deg'), frame='icrs')
        if verbose:
            print("# The center of the search: RA={:9.5f}, DEC={:9.5f}".format(ra_cen, dec_cen))
    else:
        if not isinstance(center, SkyCoord):
            raise TypeError("# The center coordinate should be a SkyCoord object")
        img_cen_ra_dec = center
        if verbose:
            print("# The center of the search: RA={:9.5f}, DEC={:9.5f}".format(
                center.ra, center.dec))

    # Width and height of the search box
    if radius is None:
        img_search_x = Quantity(pixel * (image.shape)[0] * size_buffer, u.arcsec)
        img_search_y = Quantity(pixel * (image.shape)[1] * size_buffer, u.arcsec)
        if verbose:
            print("# The width of the search: {:7.1f}".format(img_search_x))
            print("# The height of the search: {:7.1f}".format(img_search_y))
    else:
        if not isinstance(radius, Quantity):
            raise TypeError("# Searching radius needs to be an Astropy Quantity.")
        if verbose:
            print("# The searching radius is: {:7.2f}".format(radius))
    # Search for stars
    if tap_url is not None:
        with suppress_stdout():
            from astroquery.gaia import TapPlus, GaiaClass
            Gaia = GaiaClass(TapPlus(url=tap_url))

            if radius is not None:
                gaia_results = Gaia.query_object_async(
                    coordinate=img_cen_ra_dec,
                    radius=radius,
                    verbose=verbose)
            else:
                gaia_results = Gaia.query_object_async(
                    coordinate=img_cen_ra_dec,
                    width=img_search_x,
                    height=img_search_y,
                    verbose=verbose)
    else:
        with suppress_stdout():
            from astroquery.gaia import Gaia
            if radius is not None:
                gaia_results = Gaia.query_object_async(
                    coordinate=img_cen_ra_dec,
                    radius=radius,
                    verbose=verbose)
            else:
                if CustomQuery is None:
                    gaia_results = Gaia.query_object_async(
                        coordinate=img_cen_ra_dec,
                        width=img_search_x,
                        height=img_search_y,
                        verbose=verbose)
                elif CustomQuery is not None:
                    query = "SELECT * from gaiadr3.gaia_source where 1 = CONTAINS( POINT ( 'ICRS', ra, dec ), \
                    BOX ( 'ICRS', " + str(img_cen_ra_dec.ra.deg) + ", " + str(img_cen_ra_dec.dec.deg) + ", \
                    " + str(img_search_x.value*pixel/3600) + ", " + str(img_search_y.value*pixel/3600) + " ) ) AND " +  CustomQuery 
                    job = Gaia.launch_job_async ( query, dump_to_file = False )
                    gaia_results = job.get_results()
                    
    if gaia_results:
        # Convert the (RA, Dec) of stars into pixel coordinate
        ra_gaia = np.asarray(gaia_results['ra'])
        dec_gaia = np.asarray(gaia_results['dec'])
        x_gaia, y_gaia = wcs.wcs_world2pix(ra_gaia, dec_gaia, 0)

        # Generate mask for each star
        rmask_gaia_arcsec = mask_a * np.exp(
            -gaia_results['phot_g_mean_mag'] / mask_b)

        # Update the catalog
        gaia_results.add_column(Column(data=x_gaia, name='x_pix'))
        gaia_results.add_column(Column(data=y_gaia, name='y_pix'))
        gaia_results.add_column(
            Column(data=rmask_gaia_arcsec, name='rmask_arcsec'))
        return gaia_results

    return None

def gaia_star_mask(img, wcs, CustomQuery, pix=0.168, mask_a=694.7, mask_b=4.04,
                   size_buffer=1.4, gaia_bright=18.0, factor_b=1.3, factor_f=1.9):
    """Find stars using Gaia and mask them out if necessary.

    Using the stars found in the GAIA TAP catalog, we build a bright star mask following
    similar procedure in Coupon et al. (2017).

    We separate the GAIA stars into bright (G <= 18.0) and faint (G > 18.0) groups, and
    apply different parameters to build the mask.
    """        
    gaia_stars = image_gaia_stars(img, wcs, CustomQuery,pixel=pix, mask_a=mask_a, mask_b=mask_b, verbose=False,
        size_buffer=size_buffer)

    # Make a mask image
    msk_star = np.zeros(img.shape).astype('uint8')

    if gaia_stars is not None:
        gaia_b = gaia_stars[gaia_stars['phot_g_mean_mag'] <= gaia_bright]
        sep.mask_ellipse(msk_star, gaia_b['x_pix'], gaia_b['y_pix'],
                         gaia_b['rmask_arcsec'] / factor_b / pix,
                         gaia_b['rmask_arcsec'] / factor_b / pix, 0.0, r=1.0)

        gaia_f = gaia_stars[gaia_stars['phot_g_mean_mag'] > gaia_bright]
        sep.mask_ellipse(msk_star, gaia_f['x_pix'], gaia_f['y_pix'],
                         gaia_f['rmask_arcsec'] / factor_f / pix,
                         gaia_f['rmask_arcsec'] / factor_f / pix, 0.0, r=1.0)

        return gaia_stars, msk_star
    
def gaia2mask(filename, scale = 0.55, CustomQuery = None, return_GaiaResults = False, write_MaskImage = True):
    # read and scale the image
    ext = os.path.splitext(filename)[1]
    if ext == '.fz':
        hdu = fits.open(filename)[1]
    else:
        hdu = fits.open(filename)[0]    

    data = hdu.data
    header = hdu.header
    wcs = WCS(header)
    
    gaia_stars, msk_star = gaia_star_mask ( data, wcs, pix = scale, mask_a = 694.7, mask_b = 4, 
                                            size_buffer = 1.8, gaia_bright = 24.0, factor_b = 1.4, factor_f = 1.9,
                                          CustomQuery = CustomQuery)

    if write_MaskImage:
        mask = np.multiply(data, (~msk_star.astype(bool)))   
        mask_file = os.path.splitext(filename)[0] + '_gaia2mask.fits'
        fits.writeto ( mask_file, mask, hdu.header, overwrite = True ) 
        print ( '[+] File created: ', mask_file )
    
    if return_GaiaResults:
        return gaia_stars