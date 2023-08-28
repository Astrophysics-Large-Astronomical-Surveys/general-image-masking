#!/usr/bin/env python
# encoding: utf-8
"""Database query."""

import os
import sys
from contextlib import contextmanager

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse as mpl_ellip

from astropy import units as u
from astropy.table import Column
from astropy.units import Quantity
from astropy.coordinates import SkyCoord

from kungpao.display import display_single

plt.rc('text', usetex=True)

__all__ = ['image_gaia_stars', 'radec_gaia_stars']


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


def image_gaia_stars(image, wcs, radius=None, center=None, pixel=0.168,
                     mask_a=694.7, mask_b=4.04,
                     verbose=False, visual=False, size_buffer=1.4,
                     tap_url=None, img_size=(8, 8)):
    """Search for bright stars using GAIA catalog.

    TODO:
        Should be absorbed by the object for image later.

    TODO:
        Should have a version that just uses the local catalog.
    """
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
                gaia_results = Gaia.query_object_async(
                    coordinate=img_cen_ra_dec,
                    width=img_search_x,
                    height=img_search_y,
                    verbose=verbose)

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

        if visual:
            fig = plt.figure(figsize=img_size)
            ax1 = fig.add_subplot(111)

            ax1 = display_single(image, ax=ax1)
            # Plot an ellipse for each object
            for star in gaia_results:
                smask = mpl_ellip(
                    xy=(star['x_pix'], star['y_pix']),
                    width=(2.0 * star['rmask_arcsec'] / pixel),
                    height=(2.0 * star['rmask_arcsec'] / pixel),
                    angle=0.0)
                smask.set_facecolor('coral')
                smask.set_edgecolor('coral')
                smask.set_alpha(0.3)
                ax1.add_artist(smask)

            # Show stars
            ax1.scatter(
                gaia_results['x_pix'],
                gaia_results['y_pix'],
                c='orangered',
                s=100,
                alpha=0.9,
                marker='+')

            ax1.set_xlim(0, image.shape[0])
            ax1.set_ylim(0, image.shape[1])

            return gaia_results, fig

        return gaia_results

    return None

def radec_gaia_stars(ra, dec, radius=1.0, width=None, height=None, verbose=False,
                     tap_url=None):
    """Search for bright stars using GAIA catalog.

    TODO:
        Should be absorbed by the object for image later.

    TODO:
        Should have a version that just uses the local catalog.
    """
    # Central coordinate
    radec = SkyCoord(
        ra, dec, unit=('deg', 'deg'), frame='icrs')

    # Default searching radius is 1.0 deg
    r_search = Quantity(radius, u.degree)

    if width is not None:
        w_search = Quantity(width, u.degree)
        if height is None:
            # Search in a square
            h_search = w_search
        else:
            h_search = Quantity(height, u.degree)

    # Search for stars
    if tap_url is not None:
        from astroquery.gaia import TapPlus, GaiaClass
        Gaia = GaiaClass(TapPlus(url=tap_url))
    else:
        from astroquery.gaia import Gaia

    if width is not None:
        gaia_results = Gaia.query_object_async(
            coordinate=radec,
            width=w_search,
            height=h_search,
            verbose=verbose)
    else:
        gaia_results = Gaia.query_object_async(
            coordinate=radec,
            radius=r_search,
            verbose=verbose)

    return gaia_results