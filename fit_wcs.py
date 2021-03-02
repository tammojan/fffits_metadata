#!/usr/bin/env python

# Copyright (C) 2021 Cees Bassa
# SPDX-License-Identifier: GPL-3.0-or-later
"""Fit WCS to measurements"""

import numpy as np

from astropy import wcs

import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def solve_linear_equation(a, b):
    q, r = np.linalg.qr(a)
    y = np.dot(q.T, b)
    x = np.linalg.solve(r, y)
    return x


def fit_wcs(x, y, ra, dec, x0, y0, ra0, dec0, order, projection="TAN"):
    """
    Fit an astropy WCS object to given coordinates

    Args:
      x (array[float]): Pixel x coordinates
      y (array[float]): Pixel y coordinates
      ra (array[float]: RA coordinates in degrees
      dec (array[float]: Dec coordinates in degrees
      x0 (float): center pixel x
      y0 (float): center pixel y
      ra0 (float): center pixel RA in degrees
      dec0 (float) center pixel DEC in degrees
      order: order of imaging polynomials

    Returns:
      WCS instance
    """
    dx, dy = x - x0, y - y0
    ixs, iys = np.meshgrid(np.arange(order + 1), np.arange(order + 1))
    c = ixs + iys <= order
    ix, iy = ixs[c], iys[c]
    a = np.array([dx**ix[i] * dy**iy[i] for i in range(len(iy))]).T

    for k in range(5):
        w = wcs.WCS(naxis=2)
        w.wcs.ctype = ["RA---" + projection, "DEC--" + projection]
        w.wcs.cd = [[1.0, 0.0], [0.0, 1.0]]
        w.wcs.crval = [ra0, dec0]
        w.wcs.crpix = [0.0, 0.0]
        rx, ry = w.wcs_world2pix(np.stack((ra, dec), axis=-1), 1).T
        ax = solve_linear_equation(a, rx)
        ay = solve_linear_equation(a, ry)
        ra0, dec0 = w.wcs_pix2world(([[ax[0], ay[0]]]), 1)[0]
    cd = np.array([[ax[1], ax[order + 1]], [ay[1], ay[order + 1]]])
    cdinv = np.linalg.inv(cd)
    axm = np.zeros_like(ixs).astype("float32")
    aym = np.zeros_like(ixs).astype("float32")
    for i in range(len(ix)):
        if ix[i] + iy[i] >= 2:
            p = np.matmul(cdinv, np.array([ax[i], ay[i]]))
            axm[iy[i], ix[i]] = p[0]
            aym[iy[i], ix[i]] = p[1]
    w = wcs.WCS(naxis=2)
    w.wcs.ctype = ["RA---ZEA-SIP", "DEC--ZEA-SIP"]
    w.wcs.cd = cd
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [x0, y0]
    w.sip = wcs.Sip(axm.T, aym.T, None, None, w.wcs.crpix)

    return w
