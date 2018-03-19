from __future__ import (absolute_import, print_function, division)

import glob
import numpy as np

from astropy.table import Table, Column

from rebin_ndarray import bin_ndarray


def trunc_rebin(x, rfac):
    npts = len(x)
    # truncate the the array so it is an integer multiple of rebinned size
    x = x[0:int(npts/rfac)*rfac]
    return bin_ndarray(x, (int(len(x)/rfac), ), operation='mean')


def save_rebin_file(in_file, out_file, rebin_fac):
    ctable = Table.read(in_file)

    # check wavelength grid
    x = ctable['WAVELENGTH']*1e-4
    indxs, = np.where((x > 0.6) & (x < 29.))
    x = x[indxs]
    delta = x[1:] - x[0:-1]
    # print(x)
    # print(x[1:]/delta)
    print(np.mean(x[1:]/delta))
    # exit()

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

    astarnames = ['1802271', '1812095', 'bd60d1753']
    gstarnames = ['p330e', 'p177d', 'snap2']
    wdstarnames = ['g191b2b', 'gd71', 'gd153']

    # rebin to a resolution of 3000
    # a and g star models at R=300,000 and wavelength grid at R=300,000
    rfac = 50
    rb_table = Table()
    for cname in np.concatenate((astarnames, gstarnames)):
        print(cname)
        cfile = glob.glob("data/%s_mod_0??.fits" % cname)
        rb_filename = cfile[0].replace('.fits', '_r3000.fits')
        save_rebin_file(cfile[0], rb_filename, rfac)
    # wd star models at R=30,000 and wavelength grid of R=60,000
    rfac = 10
    rb_table = Table()
    for cname in wdstarnames:
        print(cname)
        rb_table = Table()
        cfile = glob.glob("data/%s_mod_0??.fits" % cname)
        rb_filename = cfile[0].replace('.fits', '_r3000.fits')
        save_rebin_file(cfile[0], rb_filename, rfac)
