from __future__ import (absolute_import, print_function, division)

import glob
import numpy as np

from astropy.table import Table, Column

from rebin_ndarray import bin_ndarray
from jwst_starnames import (astarnames_prime, astarnames_second,
                            gstarnames_prime, gstarnames_second,
                            wdstarnames_prime, wdstarnames_second)


def trunc_rebin(x, rfac):
    npts = len(x)
    # truncate the the array so it is an integer multiple of rebinned size
    x = x[0:int(npts/rfac)*rfac]
    return bin_ndarray(x, (int(len(x)/rfac), ), operation='mean')


def get_spec_waveres(ctable):
    # check wavelength grid
    x = ctable['WAVELENGTH']*1e-4
    indxs, = np.where((x > 0.6) & (x < 29.))
    x = x[indxs]
    delta = x[1:] - x[0:-1]
    return np.mean(x[1:]/delta)


def save_rebin_file(in_file, out_file, target_waveres):
    ctable = Table.read(in_file)

    actual_waveres = get_spec_waveres(ctable)
    rebin_fac = int(actual_waveres/target_waveres)
    print('actual: ', rebin_fac, actual_waveres)

    rb_table = Table()
    rb_table['WAVELENGTH'] = Column(trunc_rebin(ctable['WAVELENGTH']*1e-4,
                                                rebin_fac),
                                    description='wavelength [microns]',
                                    unit='micron')
    # should add that the units are flam, but the writer does not
    # recognize this unit - probably some astropy units magic needed
    rb_table['FLUX'] = Column(trunc_rebin(ctable['FLUX'], rebin_fac),
                              description='flux')
    if 'CONTINUUM' in rb_table.colnames:
        rb_table['CONTINUUM'] = Column(trunc_rebin(ctable['CONTINUUM'],
                                                   rebin_fac),
                                       description='continuum')
    rb_table.write(out_file, overwrite=True)


if __name__ == '__main__':

    astarnames = astarnames_prime + astarnames_second
    gstarnames = gstarnames_prime + gstarnames_second
    wdstarnames = wdstarnames_prime + wdstarnames_second

    # rebin to a resolution of 3000, wavelength grid res of 6000
    rfac = 50
    for cname in np.concatenate((astarnames, gstarnames, wdstarnames)):
        print(cname)
        cfile = glob.glob("data/%s_mod_0??.fits" % cname)
        rb_filename = cfile[0].replace('.fits', '_r3000.fits')
        save_rebin_file(cfile[0], rb_filename, 6000.0)
