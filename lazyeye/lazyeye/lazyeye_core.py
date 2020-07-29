
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

import astroquery as aq
from astroquery.sdss import SDSS
from astroquery.cadc import Cadc

from astropy import units as un
from astropy.table import vstack, join, Table
from astropy.io import fits
import astropy.coordinates as acoords


__all__ = ['sdss_find', 'sdss_images']

cadc = Cadc()

def sdss_find(co, radi = 2, drn = 16):
        """
        sdss_find takes astropy coordinates, a buffer radius and a data release to search and returns basic information on the object to be passed to other queries.

        This function makes use of and is similar in usage case to astroquery.sdss.SDSS.query_region() but utilizes preset return fields for further search.

        Parameters
        ----------
        co : str or `astropy.coordinates` object or list of
            coordinates or `~astropy.table.Column` of coordinates The
            target(s) around which to search. It may be specified as a
            string in which case it is resolved using online services or as
            the appropriate `astropy.coordinates` object. ICRS coordinates
            may also be entered as strings as specified in the
            `astropy.coordinates` module.

            Example:
            ra = np.array([110.00, 120.00, 123.00])
            dec = np.array([0.11, 0.12, 0.13])
            co = SkyCoord(ra, dec, frame='icrs', unit='deg')
        radius : str or `~astropy.units.Quantity` object, optional The
            string must be parsable by `~astropy.coordinates.Angle`. The
            appropriate `~astropy.units.Quantity` object from
            `astropy.units` may also be used. Defaults to 2 arcsec.
        data_release : int
            The data release of the SDSS to use.

        Examples
        --------
        >>> from astroquery.sdss import SDSS
        >>> from astropy import coordinates as coords
        >>> co = coords.SkyCoord('0h8m05.63s +14d50m23.3s')
        >>> result = lazyeye.sdss_find(co, 3, 15)
        >>> print(result[:5])
              ra           dec             objid        run  rerun camcol field
        ------------- ------------- ------------------- ---- ----- ------ -----
        2.02344282607 14.8398204075 1237653651835781245 1904   301      3   163
        2.02344283666 14.8398204143 1237653651835781244 1904   301      3   163
        2.02344596595 14.8398237229 1237652943176138867 1739   301      3   315
        2.02344596303 14.8398237521 1237652943176138868 1739   301      3   315
        2.02344772021 14.8398201105 1237653651835781243 1904   301      3   163

        Returns
        -------
        result : `~astropy.table.Table`
            The result of the query as a `~astropy.table.Table` object.

        """
        spec_fields = ['specObjID', 'z', 'class', 'subClass'] 

        phot_fields = ['ra', 'dec', 'ObjID', 'type', 'b', 'l', 'specObjID']

        result = SDSS.query_region(co, spectro = True, specobj_fields = spec_fields, photoobj_fields = phot_fields, data_release = drn, radius = radi * un.arcsec)

        return result

def sdss_images(co, band = 'r'):
        """
        sdss_image

        Parameters
        ----------
        co : str or `astropy.coordinates` object or list of
            coordinates or `~astropy.table.Column` of coordinates The

        band: str, list
            Could be individual band, or list of bands. Options: 'u', 'g', 'r', 'i', or 'z'.

        Examples
        --------

        Returns
        -------
        result : `~astropy.table.Table`
            The result of the query as a `~astropy.table.Table` object.

        """
        
        # SDSS.query_region(co).show_in_browser()
        images = SDSS.get_images(matches = SDSS.query_region(co), band = band)

        return images

def mega_images(co):
    images = Cadc.get_images(co, radius = 5 * un.arcsec, collection=None, show_progress=True)

    return images