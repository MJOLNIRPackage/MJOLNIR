#!/usr/bin/env python
#
# resolution ellipsoid calculations
# 
# @author Tobias Weber <tweber@ill.fr>
# @date mar-2019
# @license GPLv2
#
# @desc for algorithm: [eck14] G. Eckold and O. Sobolev, NIM A 752, pp. 54-64 (2014)
# @desc see covariance calculations: https://code.ill.fr/scientific-software/takin/mag-core/blob/master/tools/tascalc/cov.py
#

# requires numpy version >= 1.10
import numpy as np
import numpy.linalg as la
import os

#np.set_printoptions(floatmode = "fixed",  precision = 4)


sig2fwhm = 2.*np.sqrt(2.*np.log(2.))


#
# volume of the ellipsoid
#
def ellipsoid_volume(mat):
    det = np.abs(la.det(mat))

    return 4./3. * np.pi * np.sqrt(1./det)


#
# projects along one axis of the quadric
# (see [eck14], equ. 57)
#
def quadric_proj(_E, idx):
    E = np.delete(np.delete(_E, idx, axis=0), idx, axis=1)
    if np.abs(_E[idx, idx]) < 1e-8:
        return E

    v = 0.5 * (_E[idx,:] + _E[:,idx])
    vv = np.outer(v, v) / _E[idx, idx]
    vv = np.delete(np.delete(vv, idx, axis=0), idx, axis=1)

    return E - vv


#
# projects linear part of the quadric
# (see [eck14], equ. 57)
#
def quadric_proj_vec(vec, _E, idx):
    _col = _E[:,idx]
    col = np.delete(_col, idx, axis=0)
    if np.abs(_col[idx]) < 1e-8:
        return col

    v = np.delete(vec, idx, axis=0)
    v = v - col*vec[idx]/_col[idx]

    return v


#
# coherent fwhm widths
#
def calc_coh_fwhms(reso):
    vecFwhms = []
    for i in range(len(reso)):
        vecFwhms.append(sig2fwhm / np.sqrt(reso[i,i]))

    return np.array(vecFwhms)


#
# incoherent fwhm width
#
def calc_incoh_fwhm(reso):
	Qres_proj = quadric_proj(reso, 2)
	Qres_proj = quadric_proj(Qres_proj, 1)
	Qres_proj = quadric_proj(Qres_proj, 0)

	return 1./np.sqrt(np.abs(Qres_proj[0,0])) * sig2fwhm


#
# calculates the characteristics of a given ellipse by principal axis trafo
#
def descr_ellipse(quadric):
	[ evals, evecs ] = la.eig(quadric)

	fwhms = 1./np.sqrt(np.abs(evals)) * sig2fwhm

	angles = np.array([])
	if len(quadric) == 2:
		angles = np.array([ np.arctan2(evecs[1][0], evecs[0][0]) ])

	return [fwhms, angles/np.pi*180., evecs]


#
# describes the ellipsoid by a principal axis trafo and by 2d cuts
#
def calc_ellipses(Qres_Q, verbose=True):
	# 4d ellipsoid
	[fwhms, angles, rot] = descr_ellipse(Qres_Q)

	if verbose:
		print()
		print("Coherent-elastic fwhms: %s" % (calc_coh_fwhms(Qres_Q)))
		print("Incoherent-elastic fwhm: %.4f meV" % (calc_incoh_fwhm(Qres_Q)))
		print("Principal axes fwhms: %s" % fwhms)


	# 2d sliced ellipses
	if verbose:
		print()
	Qres_QxE = np.delete(np.delete(Qres_Q, 2, axis=0), 2, axis=1)
	Qres_QxE = np.delete(np.delete(Qres_QxE, 1, axis=0), 1, axis=1)
	[fwhms_QxE, angles_QxE, rot_QxE] = descr_ellipse(Qres_QxE)
	if verbose:
		print("2d Qx,E slice fwhms and slope angle: %s, %.4f" % (fwhms_QxE, angles_QxE[0]))

	Qres_QyE = np.delete(np.delete(Qres_Q, 2, axis=0), 2, axis=1)
	Qres_QyE = np.delete(np.delete(Qres_QyE, 0, axis=0), 0, axis=1)
	[fwhms_QyE, angles_QyE, rot_QyE] = descr_ellipse(Qres_QyE)
	if verbose:
		print("2d Qy,E slice fwhms and slope angle: %s, %.4f" % (fwhms_QyE, angles_QyE[0]))

	Qres_QzE = np.delete(np.delete(Qres_Q, 1, axis=0), 1, axis=1)
	Qres_QzE = np.delete(np.delete(Qres_QzE, 0, axis=0), 0, axis=1)
	[fwhms_QzE, angles_QzE, rot_QzE] = descr_ellipse(Qres_QzE)
	if verbose:
		print("2d Qz,E slice fwhms and slope angle: %s, %.4f" % (fwhms_QzE, angles_QzE[0]))

	Qres_QxQy = np.delete(np.delete(Qres_Q, 3, axis=0), 3, axis=1)
	Qres_QxQy = np.delete(np.delete(Qres_QxQy, 2, axis=0), 2, axis=1)
	[fwhms_QxQy, angles_QxQy, rot_QxQy] = descr_ellipse(Qres_QxQy)
	if verbose:
		print("2d Qx,Qy slice fwhms and slope angle: %s, %.4f" % (fwhms_QxQy, angles_QxQy[0]))


	# 2d projected ellipses
	if verbose:
		print()
	Qres_QxE_proj = np.delete(np.delete(Qres_Q, 2, axis=0), 2, axis=1)
	Qres_QxE_proj = quadric_proj(Qres_QxE_proj, 1)
	[fwhms_QxE_proj, angles_QxE_proj, rot_QxE_proj] = descr_ellipse(Qres_QxE_proj)
	if verbose:
		print("2d Qx,E projection fwhms and slope angle: %s, %.4f" % (fwhms_QxE_proj, angles_QxE_proj[0]))

	Qres_QyE_proj = np.delete(np.delete(Qres_Q, 2, axis=0), 2, axis=1)
	Qres_QyE_proj = quadric_proj(Qres_QyE_proj, 0)
	[fwhms_QyE_proj, angles_QyE_proj, rot_QyE_proj] = descr_ellipse(Qres_QyE_proj)
	if verbose:
		print("2d Qy,E projection fwhms and slope angle: %s, %.4f" % (fwhms_QyE_proj, angles_QyE_proj[0]))

	Qres_QzE_proj = np.delete(np.delete(Qres_Q, 1, axis=0), 1, axis=1)
	Qres_QzE_proj = quadric_proj(Qres_QzE_proj, 0)
	[fwhms_QzE_proj, angles_QzE_proj, rot_QzE_proj] = descr_ellipse(Qres_QzE_proj)
	if verbose:
		print("2d Qz,E projection fwhms and slope angle: %s, %.4f" % (fwhms_QzE_proj, angles_QzE_proj[0]))

	Qres_QxQy_proj = quadric_proj(Qres_Q, 3)
	Qres_QxQy_proj = np.delete(np.delete(Qres_QxQy_proj, 2, axis=0), 2, axis=1)
	[fwhms_QxQy_proj, angles_QxQy_proj, rot_QxQy_proj] = descr_ellipse(Qres_QxQy_proj)
	if verbose:
		print("2d Qx,Qy projection fwhms and slope angle: %s, %.4f" % (fwhms_QxQy_proj, angles_QxQy_proj[0]))


	results = {
		"fwhms_QxE" : fwhms_QxE, "rot_QxE" : rot_QxE, 
		"fwhms_QyE" : fwhms_QyE, "rot_QyE" : rot_QyE,
		"fwhms_QzE" : fwhms_QzE, "rot_QzE" : rot_QzE, 
		"fwhms_QxQy" : fwhms_QxQy,  "rot_QxQy" : rot_QxQy,
		"fwhms_QxE_proj" : fwhms_QxE_proj,  "rot_QxE_proj" : rot_QxE_proj, 
		"fwhms_QyE_proj" : fwhms_QyE_proj, "rot_QyE_proj" : rot_QyE_proj,
		"fwhms_QzE_proj" : fwhms_QzE_proj, "rot_QzE_proj" : rot_QzE_proj, 
		"fwhms_QxQy_proj" : fwhms_QxQy_proj, "rot_QxQy_proj" : rot_QxQy_proj,
	}

	return results



#
# shows the 2d ellipses
#
def plot_ellipses(ellis, verbose=True, plot_results=True, file="", dpi=600, ellipse_points=128):
	import mpl_toolkits.mplot3d as mplot3d
	import matplotlib.pyplot as plot

	ellfkt = lambda rad, rot, phi : \
		np.dot(rot, np.array([ rad[0]*np.cos(phi), rad[1]*np.sin(phi) ]))


	phi = np.linspace(0, 2.*np.pi, ellipse_points)

	ell_QxE = ellfkt(ellis["fwhms_QxE"]*0.5, ellis["rot_QxE"], phi)
	ell_QyE = ellfkt(ellis["fwhms_QyE"]*0.5, ellis["rot_QyE"], phi)
	ell_QzE = ellfkt(ellis["fwhms_QzE"]*0.5, ellis["rot_QzE"], phi)
	ell_QxQy = ellfkt(ellis["fwhms_QxQy"]*0.5, ellis["rot_QxQy"], phi)

	ell_QxE_proj = ellfkt(ellis["fwhms_QxE_proj"]*0.5, ellis["rot_QxE_proj"], phi)
	ell_QyE_proj = ellfkt(ellis["fwhms_QyE_proj"]*0.5, ellis["rot_QyE_proj"], phi)
	ell_QzE_proj = ellfkt(ellis["fwhms_QzE_proj"]*0.5, ellis["rot_QzE_proj"], phi)
	ell_QxQy_proj = ellfkt(ellis["fwhms_QxQy_proj"]*0.5, ellis["rot_QxQy_proj"], phi)


	# Qpara, E axis
	fig = plot.figure()
	subplot_QxE = fig.add_subplot(221)
	subplot_QxE.set_xlabel("Qpara (1/A)")
	subplot_QxE.set_ylabel("E (meV)")
	subplot_QxE.plot(ell_QxE[0], ell_QxE[1], c="black", linestyle="dashed")
	subplot_QxE.plot(ell_QxE_proj[0], ell_QxE_proj[1], c="black", linestyle="solid")

	# Qperp, E axis
	subplot_QyE = fig.add_subplot(222)
	subplot_QyE.set_xlabel("Qperp (1/A)")
	subplot_QyE.set_ylabel("E (meV)")
	subplot_QyE.plot(ell_QyE[0], ell_QyE[1], c="black", linestyle="dashed")
	subplot_QyE.plot(ell_QyE_proj[0], ell_QyE_proj[1], c="black", linestyle="solid")

	# Qup, E axis
	subplot_QzE = fig.add_subplot(223)
	subplot_QzE.set_xlabel("Qup (1/A)")
	subplot_QzE.set_ylabel("E (meV)")
	subplot_QzE.plot(ell_QzE[0], ell_QzE[1], c="black", linestyle="dashed")
	subplot_QzE.plot(ell_QzE_proj[0], ell_QzE_proj[1], c="black", linestyle="solid")

	# Qpara, Qperp axis
	subplot_QxQy = fig.add_subplot(224)
	subplot_QxQy.set_xlabel("Qpara (1/A)")
	subplot_QxQy.set_ylabel("Qperp (1/A)")
	subplot_QxQy.plot(ell_QxQy[0], ell_QxQy[1], c="black", linestyle="dashed")
	subplot_QxQy.plot(ell_QxQy_proj[0], ell_QxQy_proj[1], c="black", linestyle="solid")
	plot.tight_layout()


	# 3d plot
	fig3d = plot.figure()
	subplot3d = fig3d.add_subplot(111, projection="3d")

	subplot3d.set_xlabel("Qpara (1/A)")
	subplot3d.set_ylabel("Qperp (1/A)")
	subplot3d.set_zlabel("E (meV)")

	# xE
	subplot3d.plot(ell_QxE[0], ell_QxE[1], zs=0., zdir="y", c="black", linestyle="dashed")
	subplot3d.plot(ell_QxE_proj[0], ell_QxE_proj[1], zs=0., zdir="y", c="black", linestyle="solid")
	# yE
	subplot3d.plot(ell_QyE[0], ell_QyE[1], zs=0., zdir="x", c="black", linestyle="dashed")
	subplot3d.plot(ell_QyE_proj[0], ell_QyE_proj[1], zs=0., zdir="x", c="black", linestyle="solid")
	# xy
	subplot3d.plot(ell_QxQy[0], ell_QxQy[1], zs=0., zdir="z", c="black", linestyle="dashed")
	subplot3d.plot(ell_QxQy_proj[0], ell_QxQy_proj[1], zs=0., zdir="z", c="black", linestyle="solid")


	# save plots to files
	if file != "":
		splitext = os.path.splitext(file)
		file3d = splitext[0] + "_3d" + splitext[1]

		if verbose:
			print("Saving 2d plot to \"%s\"." % file)
			print("Saving 3d plot to \"%s\"." % file3d)
		fig.savefig(file, dpi=dpi)
		fig3d.savefig(file3d, dpi=dpi)


	# show plots
	if plot_results:
		plot.show()

	return fig,fig3d