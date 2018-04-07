from __future__ import (absolute_import, print_function, division)

import glob

from astropy.table import Table


def mock_spectra(starnames, starkmags, template_file):
    for cname, kmag in zip(starnames, starkmags):
        print(cname)
        ctable = Table.read(template_file)
        mfac = 10**(-0.4*kmag+0.4*template_star_kmag)
        ctable['FLUX'] *= mfac
        out_file = 'data/%s_mod_000.fits' % cname
        ctable.write(out_file, overwrite=True)


if __name__ == '__main__':

    # A stars
    astarnames = ['hd166205', 'hd042525', 'hd128998',
                  'hd101452', 'hd002811', 'hd055677']
    astarkmags = [4.26, 5.751, 5.756,
                  6.824, 7.043, 9.156]
    template_star = 'bd60d1753'
    template_star_kmag = 9.64
    template_file = glob.glob("data/%s_mod_0??.fits" % template_star)
    mock_spectra(astarnames, astarkmags, template_file[0])

    # G stars
    gstarnames = ['hd146233', 'hd186427', 'hd142331',
                  'hd167060', 'hd115169', 'ngc2420',
                  'ngc2506', 'ngc6811',
                  'gstark12.5', 'gstark13.0', 'gstark13.5', 'gstark14.0']
    gstarkmags = [3.99, 4.664, 7.125,
                  7.431, 7.707, 15.5,
                  16.1, 13.9,
                  12.5, 13.0, 13.5, 14.0]
    template_star = 'hd159222'
    template_star_kmag = 5.05
    template_file = glob.glob("data/%s_mod_0??.fits" % template_star)
    mock_spectra(gstarnames, gstarkmags, template_file[0])
