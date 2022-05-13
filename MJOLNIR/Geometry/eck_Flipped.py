#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# implementation of the eckold-sobolev algo
#
# @author Tobias Weber <tweber@ill.fr>
# @date feb-2015, oct-2019
# @license GPLv2
#
# Link to original code for TARKIN: https://doi.org/10.5281/zenodo.4117437
#
# @desc for algorithm: [eck14] G. Eckold and O. Sobolev, NIM A 752, pp. 54-64 (2014)
# @desc for alternate R0 normalisation: [mit84] P. W. Mitchell, R. A. Cowley and S. A. Higgins, Acta Cryst. Sec A, 40(2), 152-160 (1984)
#

# requires numpy version >= 1.10
import numpy as np
import numpy.linalg as la
from MJOLNIR.Geometry import reso


#--------------------------------------------------------------------------
# scattering triangle
# see: https://code.ill.fr/scientific-software/takin/mag-core/blob/master/tools/tascalc/tascalc.py
# 

ksq2E = 2.072124836832


def k2lam(k):
    return 2.*np.pi / k


def get_E(ki, kf):
        return ksq2E * (ki**2. - kf**2.)


def get_scattering_angle(ki, kf, Q):
    c = (ki**2. + kf**2. - Q**2.) / (2.*ki*kf)
    return np.arccos(c)


def get_angle_ki_Q(ki, kf, Q):
    c = (ki**2. + Q**2. - kf**2.) / (2.*ki*Q)
    return np.arccos(c)


def get_angle_kf_Q(ki, kf, Q):
    c = (ki**2. - Q**2. - kf**2.) / (2.*kf*Q)
    return np.arccos(c)


def get_mono_angle(k, d):
    s = np.pi/(d*k)
    angle = np.arcsin(s)
    return angle
#--------------------------------------------------------------------------



#--------------------------------------------------------------------------
# helpers

#
# z rotation matrix
#
def rotation_matrix_3d_z(angle):
    s = np.sin(angle)
    c = np.cos(angle)
    

    return np.array([
        [c, -s, 0],
        [s,  c, 0],
        [0,  0, 1]])


def mirror_matrix(iSize, iComp):
    mat = np.identity(iSize)
    mat[iComp, iComp] = -1.

    return mat;


#
# thin lens equation: 1/f = 1/lenB + 1/lenA
#
def focal_len(lenBefore, lenAfter):
    f_inv = 1./lenBefore + 1./lenAfter
    return 1. / f_inv


#
# optimal mono/ana curvature, 
# see e.g. 
# 	- (Shirane 2002) p. 66
# 	- or nicos/nicos-core.git/tree/nicos/devices/tas/mono.py in nicos
#  - or Monochromator_curved.comp in McStas
#
def foc_curv(lenBefore, lenAfter, tt, bVert):
    f = focal_len(lenBefore, lenAfter)
    s = np.abs(np.sin(0.5*tt))

    if bVert:
        curv = 2. * f*s
    else:
        curv = 2. * f/s

    return curv
#--------------------------------------------------------------------------


#
# mono (and ana) resolution calculation
#
def get_mono_vals(src_w, src_h, mono_w, mono_h,
	dist_src_mono, dist_mono_sample,
	ki, thetam,
	coll_h_pre_mono, coll_h_pre_sample,
	coll_v_pre_mono, coll_v_pre_sample,
	mono_mosaic, mono_mosaic_v,
	inv_mono_curvh, inv_mono_curvv,
	pos_x , pos_y, pos_z,
	refl):

    # A matrix: formula 26 in [eck14]
    A = np.identity(3)

    A_t0 = 1. / mono_mosaic
    A_tx = inv_mono_curvh*dist_mono_sample / np.abs(np.sin(thetam))
    A_t1 = A_t0*A_tx

    A[0,0] = 0.5*reso.sig2fwhm**2. / ki**2. * np.tan(thetam)**2. * \
        ( (2./coll_h_pre_mono)**2. + (2*dist_src_mono/src_w)**2. + A_t0*A_t0 )

    A[0,1] = A[1,0] = 0.5*reso.sig2fwhm**2. / ki**2. * np.tan(thetam) * \
        ( + 2.*(1./coll_h_pre_mono)**2. + 2.*dist_src_mono*(dist_src_mono-dist_mono_sample)/src_w**2. + \
            A_t0**2. - A_t0*A_t1)

    A[1,1] = 0.5*reso.sig2fwhm**2. / ki**2. * \
    ( (1./coll_h_pre_mono)**2. + (1./coll_h_pre_sample)**2. \
        + ((dist_src_mono-dist_mono_sample)/src_w)**2. \
        + (dist_mono_sample/(mono_w*np.abs(np.sin(thetam))))**2. \
        + A_t0*A_t0 - 2.*A_t0*A_t1 + A_t1*A_t1)



    # Av matrix: formula 38 in [eck14]
    # some typos in paper leading to the (false) result of a better Qz resolution when focusing
    # => trying to match terms in Av with corresponding terms in A
    # corresponding pre-mono terms commented out in Av, as they are not considered there
    Av = np.identity(2)

    Av_t0 = 0.5 / (mono_mosaic_v*np.abs(np.sin(thetam)))
    Av_t1 = inv_mono_curvv*dist_mono_sample / mono_mosaic_v

    Av[0,0] = 0.5*reso.sig2fwhm**2. / ki**2. * \
        ( (1./coll_v_pre_sample)**2. + (dist_mono_sample/src_h)**2. + (dist_mono_sample/mono_h)**2. + \
    	Av_t0**2. - 2.*Av_t0*Av_t1 + Av_t1**2. ) 	# typo/missing in paper?

    Av[0,1] = Av[1,0] = 0.5*reso.sig2fwhm**2. / ki**2. * \
        ( dist_src_mono*dist_mono_sample/src_h**2. - Av_t0*Av_t0 + Av_t0*Av_t1 )

    Av[1,1] = 0.5*reso.sig2fwhm**2. / ki**2. * \
        ( (1./(coll_v_pre_mono))**2. + (dist_src_mono/src_h)**2. + Av_t0**2. )



    # B vector: formula 27 in [eck14]
    B = np.array([0,0,0])
    B_t0 = inv_mono_curvh / (mono_mosaic*mono_mosaic*np.abs(np.sin(thetam)))

    B[0] = reso.sig2fwhm**2. * pos_y / ki * np.tan(thetam) * \
        ( 2.*dist_src_mono / src_w**2. + B_t0 )

    B[1] = reso.sig2fwhm**2. * pos_y / ki * \
    ( - dist_mono_sample / (mono_w*np.abs(np.sin(thetam)))**2. + \
        B_t0 - B_t0 * inv_mono_curvh*dist_mono_sample / np.abs(np.sin(thetam)) + \
        (dist_src_mono-dist_mono_sample) / src_w**2. )



    # Bv vector: formula 39 in [eck14]
    Bv = np.array([0,0])

    Bv_t0 = inv_mono_curvv/mono_mosaic_v**2

    # typo in paper?
    Bv[0] = (-1.) *  reso.sig2fwhm**2. * pos_z / ki * \
        ( dist_mono_sample / mono_h**2. + dist_mono_sample / src_h**2. + \
            Bv_t0 * inv_mono_curvv*dist_mono_sample - 0.5*Bv_t0 / np.abs(np.sin(thetam)) )

    # typo in paper?
    Bv[1] = (-1.) * reso.sig2fwhm**2. * pos_z / ki * \
        ( dist_src_mono / (src_h*src_h) + 0.5*Bv_t0/np.abs(np.sin(thetam)) )




    # C scalar: formula 28 in [eck14]
    C = 0.5*reso.sig2fwhm**2. * pos_y**2. * \
	    ( 1./src_w**2. + (1./(mono_w*np.abs(np.sin(thetam))))**2. + \
		    (inv_mono_curvh/(mono_mosaic * np.abs(np.sin(thetam))))**2. )

    # Cv scalar: formula 40 in [eck14] 
    Cv = 0.5*reso.sig2fwhm**2. * pos_z**2. * \
        ( 1./src_h**2. + 1./mono_h**2. + (inv_mono_curvv/mono_mosaic_v)**2. )



    # z components, [eck14], equ. 42
    A[2,2] = Av[0,0] - Av[0,1]*Av[0,1]/Av[1,1]
    B[2] = Bv[0] - Bv[1]*Av[0,1]/Av[1,1]
    D = Cv - 0.25*Bv[1]/Av[1,1]


    # [eck14], equ. 54
    therefl = refl * np.sqrt(np.pi / Av[1,1]) # typo in paper?
    
    return [ A, B, C, D, therefl ]




#
# Eckold algorithm combining the mono and ana resolutions
#
def calc_eck(param):
    twotheta = param["twotheta"] * param["sample_sense"]
    thetam = param["thetam"] * param["mono_sense"]
    thetaa = param["thetaa"] * param["ana_sense"]
    ki_Q = param["angle_ki_Q"] * param["sample_sense"]
    kf_Q = param["angle_kf_Q"] * param["sample_sense"]

    ki = param["ki"]
    kf = param["kf"]
    E = param["E"]
    Q = param["Q"]


    # --------------------------------------------------------------------
    # mono/ana focus
    mono_curvh = param["mono_curvh"]
    mono_curvv = param["mono_curvv"]
    ana_curvh = param["ana_curvh"]
    ana_curvv = param["ana_curvv"]

    if param["mono_is_optimally_curved_h"]:
        mono_curvh = foc_curv(param["dist_src_mono"], param["dist_mono_sample"], np.abs(2.*thetam), False)
    if param["mono_is_optimally_curved_v"]: 
        mono_curvv = foc_curv(param["dist_src_mono"], param["dist_mono_sample"], np.abs(2.*thetam), True)
    if param["ana_is_optimally_curved_h"]: 
        ana_curvh = foc_curv(param["dist_sample_ana"], param["dist_ana_det"], np.abs(2.*thetaa), False)
    if param["ana_is_optimally_curved_v"]: 
        ana_curvv = foc_curv(param["dist_sample_ana"], param["dist_ana_det"], np.abs(2.*thetaa), True)

    inv_mono_curvh = 0.
    inv_mono_curvv = 0.
    inv_ana_curvh = 0.
    inv_ana_curvv = 0.

    if param["mono_is_curved_h"]:
        inv_mono_curvh = 1./mono_curvh
    if param["mono_is_curved_v"]: 
        inv_mono_curvv = 1./mono_curvv
    if param["ana_is_curved_h"]: 
        inv_ana_curvh = 1./ana_curvh
    if param["ana_is_curved_v"]: 
        inv_ana_curvv = 1./ana_curvv
    # --------------------------------------------------------------------


    lam = k2lam(ki)

    coll_h_pre_mono = param["coll_h_pre_mono"]
    coll_v_pre_mono = param["coll_v_pre_mono"]

    if param["use_guide"]:
        coll_h_pre_mono = lam*param["guide_div_h"]
        coll_v_pre_mono = lam*param["guide_div_v"]



    # dict with results
    res = {}

    res["Q_avg"] = np.array([ Q, 0., 0., E ]).flatten()

    # -------------------------------------------------------------------------

    # - if the instruments works in kf=const mode and the scans are counted for
    #   or normalised to monitor counts no ki^3 or kf^3 factor is needed.
    # - if the instrument works in ki=const mode the kf^3 factor is needed.

    # TODO
    #tupScFact = get_scatter_factors(param.flags, param.thetam, param.ki, param.thetaa, param.kf);
    tupScFact = [1., 1., 1.]

    dmono_refl = param["dmono_refl"] * tupScFact[0]
    dana_effic = param["dana_effic"] * tupScFact[1]
    dxsec = tupScFact[2]
    #if param.mono_refl_curve:
    #    dmono_refl *= (*param.mono_refl_curve)(param.ki)
    #if param.ana_effic_curve:
    #    dana_effic *= (*param.ana_effic_curve)(param.kf)


    #--------------------------------------------------------------------------
    # mono part

    [A, B, C, D, dReflM] = get_mono_vals(
        param["src_w"], param["src_h"],
        param["mono_w"], param["mono_h"],
        param["dist_src_mono"], param["dist_mono_sample"],
        ki, thetam,
        coll_h_pre_mono, param["coll_h_pre_sample"],
        coll_v_pre_mono, param["coll_v_pre_sample"],
        param["mono_mosaic"], param["mono_mosaic_v"],
        inv_mono_curvh, inv_mono_curvv,
        param["pos_x"] , param["pos_y"], param["pos_z"],
        dmono_refl)
    #--------------------------------------------------------------------------


    #--------------------------------------------------------------------------
    # ana part
    # equ 43 in [eck14]
    
    ## UPDATED -> Flipped done 19/04/2021
    #
    # y -> z
    # z -> -y
    # w_D <-> h_D
    # w_A <-> h_A
    # alpha_3 <-> alpha_3^v
    # alpha_4 <-> alpha_4^v
    #
    pos_y2 = - param["pos_x"]*np.sin(twotheta) + param["pos_y"]*np.cos(twotheta)

    [E, F, G, H, dReflA] = get_mono_vals(
        param["det_h"], param["det_w"],
        param["ana_h"], param["ana_w"],
        param["dist_ana_det"], param["dist_sample_ana"],
        kf, -thetaa,
        param["coll_v_post_ana"], param["coll_v_post_sample"],
        param["coll_h_post_ana"], param["coll_h_post_sample"],
        param["ana_mosaic_v"], param["ana_mosaic"],
        inv_ana_curvv, inv_ana_curvh,
        param["pos_x"], param["pos_z"], -pos_y2,
        dana_effic)
    #--------------------------------------------------------------------------

    flipper = np.array([[1, 0, 0],[0, 0, -1],[0, 1, 0]])
    
    E = np.dot(flipper,np.dot(E,np.transpose(flipper)))
    F = np.dot(flipper,np.dot(F,np.transpose(flipper)))

    # equ 4 & equ 53 in [eck14]
    dE = (ki**2. - kf**2.) / (2.*Q**2.)
    kipara = Q*(0.5+dE)
    kfpara = Q-kipara
    kperp = np.sqrt(np.abs(kipara**2. - ki**2.))
    kperp *= param["sample_sense"]


    ## UPDATED -> Flipped Not needed 
    #
    # trafo, equ 52 in [eck14]
    T = np.identity(6)
    T[0,3] = T[1,4] = T[2,5] = -1.
    T[3,0] = 2.*ksq2E * kipara
    T[3,3] = 2.*ksq2E * kfpara
    T[3,1] = 2.*ksq2E * kperp
    T[3,4] = -2.*ksq2E * kperp
    T[4,1] = T[5,2] = (0.5 - dE)
    T[4,4] = T[5,5] = (0.5 + dE)

    Tinv = la.inv(T)


    # equ 54 in [eck14]
    Dalph_i = rotation_matrix_3d_z(-ki_Q)
    Dalph_f = rotation_matrix_3d_z(-kf_Q)
    Arot = np.dot(np.dot(np.transpose(Dalph_i), A), Dalph_i)
    Erot = np.dot(np.dot(np.transpose(Dalph_f), E), Dalph_f)

    matAE = np.zeros((6,6))
    matAE[0:3, 0:3] = Arot
    matAE[3:6, 3:6] = Erot

    # U1 matrix
    # typo in paper in quadric trafo in equ 54 (top)?
    U1 = np.dot(np.dot(np.transpose(Tinv), matAE), Tinv)

    # V1 vector
    vecBF = np.zeros(6)
    vecBrot = np.dot(np.transpose(Dalph_i), B)
    vecFrot = np.dot(np.transpose(Dalph_f), F)
    vecBF[0:3] = vecBrot
    vecBF[3:6] = vecFrot
    V1 = np.dot(vecBF, Tinv)



    # --------------------------------------------------------------------------
    # integrate last 2 vars -> equs 57 & 58 in [eck14]

    U2 = reso.quadric_proj(U1, 5);
    U = reso.quadric_proj(U2, 4);

    V2 = reso.quadric_proj_vec(V1, U1, 5);
    V = reso.quadric_proj_vec(V2, U2, 4);

    W = (C + D + G + H) - 0.25*V1[5]/U1[5,5] - 0.25*V2[4]/U2[4,4]

    Z = dReflM*dReflA * np.sqrt(np.pi/np.abs(U1[5,5])) * np.sqrt(np.pi/np.abs(U2[4,4]))
    # --------------------------------------------------------------------------



    # quadratic part of quadric (matrix U)
    # careful: factor -0.5*... missing in U matrix compared to normal gaussian!
    res["reso"] = 2. * U;
    # linear (vector V) and constant (scalar W) part of quadric
    res["reso_v"] = V;
    res["reso_s"] = W;


    if param["sample_sense"] < 0.:
        # mirror Q_perp
        matMirror = mirror_matrix(len(res["reso"]), 1)
        res["reso"] = np.dot(np.dot(np.transpose(matMirror), res["reso"]), matMirror)
        res["reso_v"][1] = -res["reso_v"][1]


    # prefactor and volume
    res["res_vol"] = reso.ellipsoid_volume(res["reso"])

    res["r0"] = Z
    # missing volume prefactor to normalise gaussian,
    # cf. equ. 56 in [eck14] to  equ. 1 in [pop75] and equ. A.57 in [mit84]
    res["r0"] *= res["res_vol"] * np.pi * 3.

    res["r0"] *= np.exp(-W)
    res["r0"] *= dxsec

	# Bragg widths
    res["coherent_fwhms"] = reso.calc_coh_fwhms(res["reso"])
    res["ok"] = True

    if np.isnan(res["r0"]) or np.isinf(res["r0"]) or np.isnan(res["reso"].any()) or np.isinf(res["reso"].any()):
        res["ok"] = False

    return res



#
# test calculation
#
if __name__ == "__main__":
    verbose = True

    cm2A = 1e8
    min2rad = 1./ 60. / 180.*np.pi
    rad2deg = 180. / np.pi

    d_mono = 3.355
    d_ana = 3.355

    ki = 2*np.pi/4
    kf = 2*np.pi/4
    Q = 1.5
    E = get_E(ki, kf)

    sc_senses = [ 1., -1., 1.]

    params = {
        # scattering triangle
        "ki" : ki, "kf" : kf, "E" : E, "Q" : Q,

        # angles
        "twotheta" : get_scattering_angle(ki, kf, Q),
        "thetam" : get_mono_angle(ki, d_mono),
        "thetaa" : get_mono_angle(kf, d_ana),
        "angle_ki_Q" : get_angle_ki_Q(ki, kf, Q),
        "angle_kf_Q" : get_angle_kf_Q(ki, kf, Q),
    
        # scattering senses
        "mono_sense" : sc_senses[0],
        "sample_sense" : sc_senses[1],
        "ana_sense" : sc_senses[2],

        # distances
        "dist_src_mono" : 100. * cm2A,
        "dist_mono_sample" : 100. * cm2A,
        "dist_sample_ana" : 100. * cm2A,
        "dist_ana_det" : 100. * cm2A,

        # component sizes
        "src_w" : 10. * cm2A,
        "src_h" : 10. * cm2A,
        "mono_w" : 10. * cm2A,
        "mono_h" : 10. * cm2A,
        "det_w" : 10. * cm2A,
        "det_h" : 10. * cm2A,
        "ana_w" : 10. * cm2A,
        "ana_h" : 10. * cm2A,

        # focusing
        "mono_curvh" : 0.,
        "mono_curvv" : 0.,
        "ana_curvh" : 0.,
        "ana_curvv" : 0.,
        "mono_is_optimally_curved_h" : False,
        "mono_is_optimally_curved_v" : False,
        "ana_is_optimally_curved_h" : False,
        "ana_is_optimally_curved_v" : False,
        "mono_is_curved_h" : False,
        "mono_is_curved_v" : False,
        "ana_is_curved_h" : False,
        "ana_is_curved_v" : False,

        # collimation
        "coll_h_pre_mono" : 9999. *min2rad,
        "coll_v_pre_mono" : 9999. *min2rad,
        "coll_h_pre_sample" : 9999. *min2rad,
        "coll_v_pre_sample" : 9999. *min2rad,
        "coll_h_post_sample" : 9999. *min2rad,
        "coll_v_post_sample" : 9999. *min2rad,
        "coll_h_post_ana" : 9999. *min2rad,
        "coll_v_post_ana" : 9999. *min2rad,

        # guide
        "use_guide" : True,
        "guide_div_h" : 9999. *min2rad,
        "guide_div_v" : 9999. *min2rad,

        # mosaics
        "mono_mosaic" : 60. *min2rad,
        "mono_mosaic_v" : 60. *min2rad,
        "ana_mosaic" : 60. *min2rad,
        "ana_mosaic_v" : 60. *min2rad,

        # crystal reflectivities
        # TODO, so far always 1
        "dmono_refl" : 1.,
        "dana_effic" : 1.,

        # off-center scattering
        # WARNING: while this is calculated, it is not yet considered in the ellipse plots
        "pos_x" : 0. * cm2A,
        "pos_y" : 0. * cm2A,
        "pos_z" : 0. * cm2A,
    }


    # calculate resolution ellipsoid
    res = calc_eck(params)
    if not res["ok"]:
        print("RESOLUTION CALCULATION FAILED!")
        exit(-1)

    if verbose:
        print("2theta = %g, thetam = %g, thetaa = %g, ki_Q = %g, kf_Q = %g\n" % 
            (params["twotheta"]*rad2deg, params["thetam"]*rad2deg, params["thetaa"]*rad2deg, 
            params["angle_ki_Q"]*rad2deg, params["angle_kf_Q"]*rad2deg))
        print("R0 = %g, Vol = %g" % (res["r0"], res["res_vol"]))
        print("Resolution matrix:\n%s" % res["reso"])
        print("Resolution vector: %s" % res["reso_v"])
        print("Resolution scalar: %g" % res["reso_s"])


    # describe and plot ellipses
    ellipses = reso.calc_ellipses(res["reso"], verbose)
    reso.plot_ellipses(ellipses, verbose)

