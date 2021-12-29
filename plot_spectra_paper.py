from __future__ import absolute_import, print_function, division

import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.lines import Line2D

from astropy.table import Table, QTable
from astropy import units as u

from jwst_starnames import (
    astarnames_prime,
    astarnames_second,
    gstarnames_prime,
    gstarnames_second,
    wdstarnames_prime,
    wdstarnames_second,
)

FNU = u.erg / (u.cm ** 2 * u.s * u.Hz)
FLAM = u.erg / (u.cm ** 2 * u.s * u.AA)


def set_params(lw=1.5, universal_color="#262626", fontsize=16):
    """Configure some matplotlib rcParams.
    Parameters
    ----------
    lw : scalar
        Linewidth of plot and axis lines. Default is 1.5.
    universal_color : str, matplotlib named color, rgb tuple
        Color of text and axis spines. Default is #262626, off-black
    fontsize : scalar
        Font size in points. Default is 12
    """
    rc("font", size=fontsize)
    rc("lines", linewidth=lw)
    rc("patch", linewidth=lw, edgecolor="#FAFAFA")
    rc(
        "axes",
        linewidth=lw,
        edgecolor=universal_color,
        labelcolor=universal_color,
        axisbelow=True,
    )
    rc("image", origin="lower")
    rc("xtick.major", width=lw * 0.75)
    rc("xtick.minor", width=lw * 0.5)
    rc("xtick", color=universal_color)
    rc("ytick.major", width=lw * 0.75)
    rc("ytick.minor", width=lw * 0.5)
    rc("ytick", color=universal_color)
    rc("grid", linewidth=lw)
    rc(
        "legend",
        loc="best",
        numpoints=1,
        scatterpoints=1,
        handlelength=1.5,
        fontsize=fontsize,
        columnspacing=1,
        handletextpad=0.75,
    )


def initialize_parser():
    """For running from command line, initialize argparse with common args"""
    ftypes = [
        "png",
        "jpg",
        "jpeg",
        "pdf",
        "ps",
        "eps",
        "rgba",
        "svg",
        "tiff",
        "tif",
        "pgf",
        "svgz",
        "raw",
    ]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s",
        "--savefig",
        action="store",
        default=False,
        choices=ftypes,
        help="Save figure to a file",
    )
    parser.add_argument("--astars", help="Use A star models", action="store_true")
    parser.add_argument("--gstars", help="Use G star models", action="store_true")
    parser.add_argument("--wdstars", help="Use WD star models", action="store_true")
    parser.add_argument(
        "--inst",
        default="all",
        choices=["all", "NIRCAM", "NIRSPEC", "NIRISS", "FGS", "MIRI"],
        help="Instruments to plot",
    )
    return parser


def plot_spec_singleinst(cax, inst, modes="all"):

    cycle1_stars = [
        "1743045",
        "1802271",
        "1812095",
        "bd60d1753",
        "hd002811",
        "hd180609",
        "hd166205",
        "p330e",
        "p177d",
        "hd167060",
        "g191b2b",
        "gd71",
        "gd153",
    ]

    starfluxes = {}
    starwaves = {}
    for cname in allstarnames:
        cfile = glob.glob("data/%s_mod_0??.fits" % cname)
        rb_filename = cfile[0].replace(".fits", "_r3000.fits")
        ctable = Table.read(rb_filename)

        if cname in astarnames:
            col = "g-"
        elif cname in gstarnames:
            col = "m:"
        elif cname in wdstarnames:
            col = "b-."

        x = ctable["WAVELENGTH"].quantity
        (indxs,) = np.where((x > 0.6 * u.um) & (x < 29.0 * u.um))
        x = ctable["WAVELENGTH"][indxs].quantity
        flux = ctable["FLUX"][indxs].quantity * FLAM
        flux_mJy = flux.to(u.mJy, u.spectral_density(x))

        alpha = 0.25
        if cname in cycle1_stars:
            alpha = 0.75

        cax.plot(x, (x ** 2) * flux_mJy, col, label=cname, alpha=alpha)

        # save star fluxes and wavelengths
        starfluxes[cname] = flux_mJy
        starwaves[cname] = x

    # read in the min/max sensitivities for the modes
    sens_filename = "jwst_inst_sens_paperplot.dat"
    mmvals = QTable.read(
        sens_filename, format="ascii.commented_header", header_start=-1
    )

    # add in the units
    mmvals["wave"] *= u.micron
    mmvals["full_min"] *= u.mJy
    mmvals["full_max"] *= u.mJy
    mmvals["sub_max"] *= u.mJy

    # plot the min/max sensitivites
    uinst = inst
    ctype = {"NIRCAM": "k", "NIRSPEC": "k", "NIRISS": "k", "MIRI": "k", "FGS": "k"}
    ptype = {
        "NIRCAM": "solid",
        "NIRSPEC": "solid",
        "NIRISS": "solid",
        "MIRI": "solid",
        "FGS": "solid",
    }
    for cinst in uinst:
        (iindxs,) = np.where(mmvals["inst"] == cinst)
        umode = np.unique(mmvals["mmode"][iindxs])
        for cmode in umode:
            (mindxs,) = np.where(mmvals["mmode"][iindxs] == cmode)
            for k in mindxs:
                ll = iindxs[k]
                modename = mmvals["mmode"][ll]
                bandname = mmvals["band"][ll]
                bandmin = mmvals["full_min"][ll]
                if mmvals["sub_max"][ll] > 0:
                    bandmax = mmvals["sub_max"][ll]
                else:
                    bandmax = mmvals["full_max"][ll]
                cwave = mmvals["wave"][ll]
                ccolor = ctype[cinst]
                cline = ptype[cinst]
                if cinst == "NIRCAM":
                    if "WFSS" in modename:
                        cline = "solid"
                        bandname = modename
                    elif "CORON" in modename:
                        cline = "dashed"
                    elif "M" in bandname:
                        cline = "dashed"
                    elif "N" in bandname:
                        cline = "dotted"
                elif cinst == "NIRSPEC":
                    bandname = modename
                    if "IFU" in modename:
                        cline = "dashed"
                elif cinst == "NIRISS":
                    if "WFSS" in modename:
                        cline = "dashed"
                        bandname = modename
                    elif "SOSS" in modename:
                        cline = "dotted"
                        bandname = modename
                    elif "AMI" in modename:
                        cline = "dashdot"
                elif cinst == "MIRI":
                    if modename == "4QPM" or modename == "LYOT":
                        cline = "dashed"
                    elif "LRS" in modename:
                        cline = "dotted"
                        bandname = modename
                    elif "MRS" in modename:
                        cline = "dashdot"
                        bandname = modename

                if modes == "all" or modename in modes:
                    cax.plot(
                        [cwave.value, cwave.value],
                        np.array([bandmax.value, bandmin.value]) * cwave.value ** 2,
                        color=ccolor,
                        linestyle=cline,
                        linewidth=2.0,
                    )
                    cax.text(
                        cwave.value,
                        bandmax.value * cwave.value ** 2,
                        bandname,
                        rotation=45.0,
                        fontsize=9.0,
                        horizontalalignment="center",
                    )

    cax.set_xscale("linear")
    cax.set_xlim(0.6, 5.0)
    cax.set_yscale("log")
    cax.set_ylim(1e-4, 1e7)
    # cax.set_ylim(a_yrange)
    cax.set_ylabel(r"$\lambda^2 F(\nu)$", fontsize=12)
    # cax.set_ylabel(r"$\lambda^2 F(nu)$ [mJy $\mu m^2$]")
    # cax.legend()


if __name__ == "__main__":

    parser = initialize_parser()
    args = parser.parse_args()

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

    xsize = 12.0
    ysize = 14.0
    fig, ax = plt.subplots(nrows=5, ncols=1, figsize=(xsize, ysize))
    # xsize = 10.0
    # ysize = 18.0
    # fig, ax = plt.subplots(nrows=6, ncols=1, figsize=(xsize, ysize))

    set_params(lw=1.0, fontsize=16)

    # cax = ax[0, 0]
    cax = ax[3]
    plot_spec_singleinst(cax, ["NIRCAM"], modes=["IMAGE"])
    cax.set_ylim(1e-4, 1e6)
    legend_elements = [
        Line2D([0], [0], color="k", lw=2, linestyle="solid", label="WIDE"),
        Line2D([0], [0], color="k", lw=2, linestyle="dashed", label="MEDIUM"),
        Line2D([0], [0], color="k", lw=2, linestyle="dotted", label="NARROW"),
    ]
    cax.legend(
        handles=legend_elements,
        fontsize=10,
        loc="lower right",
        title="NIRCam IMAGE",
        title_fontsize=11,
    )

    # cax = ax[1, 0]
    # cax = ax[1]
    # plot_spec_singleinst(cax, ["NIRCAM"], modes=["WFSS", "CORON"])
    # legend_elements = [
    #     Line2D([0], [0], color="k", lw=2, linestyle="solid", label="WFSS"),
    #     Line2D([0], [0], color="k", lw=2, linestyle="dashed", label="CORON"),
    # ]
    # cax.legend(
    #     handles=legend_elements,
    #     fontsize=10,
    #     loc="lower right",
    #     title="NIRCam",
    #     title_fontsize=11,
    # )

    # cax = ax[2, 0]
    cax = ax[2]
    plot_spec_singleinst(cax, ["NIRISS"])
    plot_spec_singleinst(cax, ["FGS"])
    legend_elements = [
        Line2D([0], [0], color="k", lw=2, linestyle="solid", label="IMAGE"),
        Line2D([0], [0], color="k", lw=2, linestyle="dashed", label="WFSS"),
        Line2D([0], [0], color="k", lw=2, linestyle="dotted", label="SOSS"),
        Line2D([0], [0], color="k", lw=2, linestyle="dashdot", label="AMI"),
    ]
    cax.legend(
        handles=legend_elements,
        fontsize=10,
        loc="lower right",
        title="NIRISS/FGS",
        title_fontsize=11,
    )

    # cax = ax[0, 1]
    cax = ax[1]
    plot_spec_singleinst(cax, ["NIRSPEC"])
    legend_elements = [
        Line2D([0], [0], color="k", lw=2, linestyle="solid", label="FixedSlit"),
        Line2D([0], [0], color="k", lw=2, linestyle="dashed", label="IFU"),
    ]
    cax.legend(
        handles=legend_elements,
        fontsize=10,
        loc="lower right",
        title="NIRSpec",
        title_fontsize=11,
    )

    # cax = ax[1, 1]
    cax = ax[0]
    plot_spec_singleinst(cax, ["MIRI"])
    cax.set_xlim(5.0, 29.0)
    cax.set_ylim(1e-2, 1e8)
    legend_elements = [
        Line2D([0], [0], color="k", lw=2, linestyle="solid", label="IMAGE"),
        Line2D([0], [0], color="k", lw=2, linestyle="dashed", label="CORON"),
        Line2D([0], [0], color="k", lw=2, linestyle="dotted", label="LRS"),
        Line2D([0], [0], color="k", lw=2, linestyle="dashdot", label="MRS"),
    ]
    cax.legend(
        handles=legend_elements,
        fontsize=10,
        loc="lower right",
        title="MIRI",
        title_fontsize=11,
    )

    # cax = ax[2, 1]
    cax = ax[4]
    plot_spec_singleinst(cax, ["NIRCAM"], modes=["WFSS", "CORON"])
    legend_elements = [
        Line2D([0], [0], color="k", lw=2, linestyle="solid", label="WFSS"),
        Line2D([0], [0], color="k", lw=2, linestyle="dashed", label="CORON"),
    ]
    leg1 = cax.legend(
        handles=legend_elements,
        fontsize=10,
        loc="lower right",
        title="NIRCam",
        title_fontsize=11,
    )
    plt.gca().add_artist(leg1)
    legend_elements = [
        Line2D([0], [0], color="b", lw=2, alpha=0.25, label="WD stars"),
        Line2D([0], [0], color="g", lw=2, alpha=0.25, label="A stars"),
        Line2D([0], [0], color="m", lw=2, alpha=0.25, label="G stars"),
    ]
    leg2 = cax.legend(
        handles=legend_elements,
        fontsize=10,
        loc="upper left",
        title="Calibrators",
        title_fontsize=11,
    )
    plt.gca().add_artist(leg2)
    legend_elements = [
        Line2D([0], [0], color="b", lw=2, alpha=0.75, label="WD stars"),
        Line2D([0], [0], color="g", lw=2, alpha=0.75, label="A stars"),
        Line2D([0], [0], color="m", lw=2, alpha=0.75, label="G stars"),
    ]
    cax.legend(
        handles=legend_elements,
        fontsize=10,
        loc="lower left",
        title="Cycle 1",
        title_fontsize=11,
    )
    cax.set_xlabel(r"$\lambda$ [$\mu m$]", fontsize=12)

    fig.tight_layout()

    # save the plot
    basename = "jwst_abscal_spec_paper"
    if args.astars:
        basename += "_astars"
    if args.gstars:
        basename += "_gstars"
    if args.wdstars:
        basename += "_wdstars"
    if args.savefig:
        fig.savefig("%s.%s" % (basename, args.savefig))
    else:
        plt.show()
