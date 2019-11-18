from __future__ import (absolute_import, print_function, division)

import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

from astropy.table import Table
from astropy import units as u

from jwst_starnames import (astarnames_prime, astarnames_second,
                            gstarnames_prime, gstarnames_second,
                            wdstarnames_prime, wdstarnames_second)

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
    parser.add_argument('--astars', help='Use A star models',
                        action='store_true')
    parser.add_argument('--gstars', help='Use G star models',
                        action='store_true')
    parser.add_argument('--wdstars', help='Use WD star models',
                        action='store_true')
    parser.add_argument('--primeonly', action='store_true',
                        help='Only use prime designated sources')
    parser.add_argument('--secondonly', action='store_true',
                        help='Only use secondary designated sources')
    parser.add_argument('--part1', action='store_true',
                        help='Only check the part1 limited set of modes')
    parser.add_argument("-t", "--target_obs", metavar=int, default=3,
                        help="number of target observations per mode")
    return parser


def which_observable(modewave, modemin, modemax,
                     starnames, starwaves, starfluxes):
    """
    Determine which stars can be observed in this mode

    Parameters
    ----------
    modewave : float
        wavelength of mode

    modemin, modemax: floats
        min/max fluxes observable in this mode

    starnames: list of strings
        names of the stars

    starwaves, starfluxes: dict of wave/flux vectors
        models of the star fluxes

    Returns
    -------
    obsnames : list of str
        names of the stars that are observable in this mode
    """
    obsstars = []
    for cname in starnames:
        modeflux = np.interp([modewave], starwaves[cname], starfluxes[cname])
        if modemin <= modeflux <= modemax:
            obsstars.append(cname)

    return obsstars


def which_modes(mmvals,
                starnames, starwaves, starfluxes):
    """
    Determine which modes are observable by which star

    Parameters
    ----------
    mmvals : astropy Table
        mode min/max info

    starnames: list of strings
        names of the stars

    starwaves, starfluxes: dict of wave/flux vectors
        models of the star fluxes

    Returns
    -------
    starmodes, starmodes_num: tuple of dicts
        dictonary by starname giving the modes and number observable
    """
    starmodes = {}
    starmodes_num = {}
    for k in range(len(mmvals)):
        if mmvals['sub_max'][k] > 0:
            bandmax = mmvals['sub_max'][k]
        else:
            bandmax = mmvals['full_max'][k]
        obsnames = which_observable(mmvals['wave'][k], mmvals['full_min'][k],
                                    bandmax,
                                    starnames, starwaves, starfluxes)
        for cname in obsnames:
            modeid = (mmvals['inst'][k], mmvals['mmode'][k], mmvals['band'][k])
            if cname not in starmodes.keys():
                starmodes[cname] = []
                starmodes_num[cname] = 0
            starmodes[cname].append(modeid)
            starmodes_num[cname] += 1

    return (starmodes, starmodes_num)


if __name__ == '__main__':

    parser = initialize_parser()
    args = parser.parse_args()

    if args.primeonly:
        astarnames = astarnames_prime
        gstarnames = gstarnames_prime
        wdstarnames = wdstarnames_prime
    elif args.secondonly:
        astarnames = astarnames_second
        gstarnames = gstarnames_second
        wdstarnames = wdstarnames_second
    else:
        astarnames = astarnames_prime + astarnames_second
        gstarnames = gstarnames_prime + gstarnames_second
        wdstarnames = wdstarnames_prime + wdstarnames_second

    allstarnames = []
    if args.astars:
        allstarnames += astarnames
    if args.gstars:
        allstarnames += gstarnames
    if args.wdstars:
        allstarnames += wdstarnames
    if len(allstarnames) == 0:
        allstarnames = astarnames + gstarnames + wdstarnames

    target_num_obs = int(args.target_obs)

    xsize = 15.0
    ysize = 9.0
    fig, cax = plt.subplots(figsize=(xsize, ysize))

    set_params(lw=1., fontsize=16)

    starfluxes = {}
    starwaves = {}
    for cname in allstarnames:
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

    # read in the min/max sensitivities for the modes
    if args.part1:
        sens_filename = 'jwst_inst_sens_part1.dat'
    else:
        sens_filename = 'jwst_inst_sens.dat'
    mmvals = Table.read(sens_filename,
                        format='ascii.commented_header',
                        header_start=-1)

    # create a dictionary with all keys for all the modes
    #   using a tuple of (inst, mmode, band) as the key
    mo_keys = zip(mmvals['inst'], mmvals['mmode'], mmvals['band'])
    modeobserved = {}
    modeobservedstars = {}
    modedone = {}
    for cmode in mo_keys:
        modeobserved[cmode] = 0
        modeobservedstars[cmode] = []
        modedone[cmode] = False

    # for each stars determine the modes and number observable
    indxs = list(range(len(mmvals)))
    cstarnames = allstarnames.copy()
    cstarwaves = starwaves.copy()
    cstarfluxes = starfluxes.copy()
    obsstarlist = []
    while True:
        starmodes, sm_num = which_modes(mmvals[indxs],
                                        cstarnames, cstarwaves, cstarfluxes)
        # get the starname with the most observed modes
        sname = ''
        maxobs = 0
        for cname in sm_num.keys():
            if sm_num[cname] > maxobs:
                sname = cname
                maxobs = sm_num[cname]
            if cname == 'p330e':  # prioritize p330e
                sname = cname
                maxobs = 1000

        # stop if no star covers the remaining modes
        if sname == '':
            break

        # add this star to the list for observations
        obsstarlist.append(sname)

        # tabulate that all the modes for this stars have one more obs
        for cur_mokey in starmodes[sname]:
            modeobserved[cur_mokey] += 1
            modeobservedstars[cur_mokey].append(sname)
        # remove the star from the possible stars list
        # sn_k, = np.where(cstarnames == sname)
        for k in range(len(cstarnames)):
            if sname == cstarnames[k]:
                sn_k = k
        del cstarnames[sn_k]
        del cstarwaves[sname]
        del cstarfluxes[sname]

        # check the list of modes, remove a mode if it has the target number
        for cur_mokey in modeobserved:
            if (modeobserved[cur_mokey] >= target_num_obs):
                if not modedone[cur_mokey]:
                    modedone[cur_mokey] = True
                    cinst, cmmode, cband = cur_mokey
                    dindxs, = np.where((mmvals['inst'] == cinst)
                                       & (mmvals['mmode'] == cmmode)
                                       & (mmvals['band'] == cband))
                    dindxs2, = np.where(dindxs[0] == indxs)
                    del indxs[dindxs2[0]]

        # stop if no stars left
        if len(cstarnames) <= 0:
            break

    print('star list')
    print(obsstarlist)

    mo_keys = zip(mmvals['inst'], mmvals['mmode'], mmvals['band'])
    print('%8s  %8s  %8s  %2s  %s' % ('Inst', 'MMode', 'Band', '#', 'stars'))
    for ckey in mo_keys:
        print('%8s,  %8s,  %8s,  %2i' % (ckey[0], ckey[1], ckey[2],
                                         modeobserved[ckey]),
              modeobservedstars[ckey])

    # plot the min/max sensitivites
    uinst = np.unique(mmvals['inst'])
    ctype = {'NIRCAM': 'c', 'NIRSPEC': 'm',
             'NIRISS': 'y', 'MIRI': 'k', 'FGS': 'b'}
    ptype = {'NIRCAM': 'solid', 'NIRSPEC': 'dashed',
             'NIRISS': 'dotted', 'MIRI': 'dashdot', 'FGS': 'dashdot'}
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
    cax.set_ylabel(r"$\lambda^2 F(nu)$ [mJy $\mu m^2$]")
    # cax.legend()

    fig.tight_layout()

    # save the plot
    basename = 'jwst_abscal_spec'
    if args.astars:
        basename += '_astars'
    if args.gstars:
        basename += '_gstars'
    if args.wdstars:
        basename += '_wdstars'
    if args.part1:
        basename += '_part1'
    if args.savefig:
        fig.savefig('%s.%s' % (basename, args.savefig))
    else:
        plt.show()
