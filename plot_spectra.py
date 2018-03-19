from __future__ import (absolute_import, print_function, division)

import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

from astropy.table import Table
from astropy import units as u


FNU = u.erg / (u.cm**2 * u.s * u.Hz)
FLAM = u.erg / (u.cm**2 * u.s * u.AA)


def set_params(lw=1.5, universal_color='#262626', fontsize=16):
    '''Configure some matplotlib rcParams.
    Parameters
    ----------
    lw : scalar
        Linewidth of plot and axis lines. Default is 1.5.
    universal_color : str, matplotlib named color, rgb tuple
        Color of text and axis spines. Default is #262626, off-black
    fontsize : scalar
        Font size in points. Default is 12
    '''
    rc('font', size=fontsize)
    rc('lines', linewidth=lw)
    rc('patch', linewidth=lw, edgecolor='#FAFAFA')
    rc('axes', linewidth=lw, edgecolor=universal_color,
       labelcolor=universal_color,
       axisbelow=True)
    rc('image', origin='lower')
    rc('xtick.major', width=lw*0.75)
    rc('xtick.minor', width=lw*0.5)
    rc('xtick', color=universal_color)
    rc('ytick.major', width=lw*0.75)
    rc('ytick.minor', width=lw*0.5)
    rc('ytick', color=universal_color)
    rc('grid', linewidth=lw)
    rc('legend', loc='best', numpoints=1, scatterpoints=1, handlelength=1.5,
        fontsize=fontsize, columnspacing=1, handletextpad=0.75)


def initialize_parser():
    '''For running from command line, initialize argparse with common args
    '''
    ftypes = ['png', 'jpg', 'jpeg', 'pdf', 'ps', 'eps', 'rgba',
              'svg', 'tiff', 'tif', 'pgf', 'svgz', 'raw']
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--savefig', action='store',
                        default=False, choices=ftypes,
                        help='Save figure to a file')
    return parser


if __name__ == '__main__':

    parser = initialize_parser()
    args = parser.parse_args()

    astarnames = ['1802271', '1812095', 'bd60d1753']
    gstarnames = ['p330e', 'p177d', 'snap2']
    wdstarnames = ['g191b2b', 'gd71', 'gd153']

    xsize = 15.0
    ysize = 9.0
    fig, cax = plt.subplots(figsize=(xsize, ysize))

    set_params(lw=1., fontsize=16)

    starfluxes = {}
    starwaves = {}
    for cname in np.concatenate((astarnames, gstarnames, wdstarnames)):
        cfile = glob.glob("data/%s_mod_0??.fits" % cname)
        rb_filename = cfile[0].replace('.fits', '_r3000.fits')
        ctable = Table.read(rb_filename)

        if cname in astarnames:
            col = 'g'
        elif cname in gstarnames:
            col = 'r'
        elif cname in wdstarnames:
            col = 'b'

        x = ctable['WAVELENGTH'].quantity
        indxs, = np.where((x > 0.6*u.um) & (x < 29.*u.um))
        x = ctable['WAVELENGTH'][indxs].quantity
        flux = ctable['FLUX'][indxs].quantity*FLAM
        flux_mJy = flux.to(u.mJy, u.spectral_density(x))

        cax.plot(x, (x**2)*flux_mJy,
                 col+'-', label=cname)

        # save star fluxes and wavelengths
        starfluxes[cname] = flux_mJy
        starwaves[cname] = x

    # plot the min/max sensitivites
    mmvals = Table.read('jwst_inst_sens.dat',
                        format='ascii.commented_header',
                        header_start=-1)
    uinst = np.unique(mmvals['inst'])
    ctype = {'NIRCAM': 'c', 'NIRSPEC': 'm',
             'NIRISS': 'y', 'MIRI': 'k'}
    ptype = {'NIRCAM': 'solid', 'NIRSPEC': 'dashed',
             'NIRISS': 'dotted', 'MIRI': 'dashdot'}
    for cinst in uinst:
        iindxs, = np.where(mmvals['inst'] == cinst)
        umode = np.unique(mmvals['mmode'][iindxs])
        for cmode in umode:
            mindxs, = np.where(mmvals['mmode'][iindxs] == cmode)
            for k in mindxs:
                ll = iindxs[k]
                bandname = mmvals['band'][ll]
                bandmin = mmvals['full_min'][ll]
                if mmvals['sub_max'][ll] > 0:
                    bandmax = mmvals['sub_max'][ll]
                else:
                    bandmax = mmvals['full_max'][ll]
                cwave = mmvals['wave'][ll]
                cax.plot([cwave, cwave],
                         np.array([bandmax, bandmin])*cwave**2,
                         color=ctype[cinst], linestyle=ptype[cinst],
                         linewidth=2.)

    cax.set_xscale('log')
    cax.set_xlim(0.6, 29.0)
    cax.set_yscale('log')
    # cax.set_ylim(a_yrange)
    cax.set_ylabel("$\lambda^2 F(nu)$ [mJy $\mu m^2$]")
    # cax.legend()

    fig.tight_layout()

    # save the plot
    basename = 'jwst_abscal_spec_'
    if args.savefig:
        fig.savefig('%s.%s' % (basename, args.savefig))
    else:
        plt.show()
