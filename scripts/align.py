import numpy as np


def rmsd(geomA, geomB, prealign=True):
    """ Calculate root-mean-square deviation between two sets of points.

    Args:
        geomA: First set of points, (npoint x dim) array.
        geomB: Second set of points, (npoint x dim) array.
        pre_align: Align geometries before calculating RMSD.

    Returns:
        rmsd: Distance between points, scalar.
    """
    # Check array dimensions
    workA = check_geom_array_shape(geomA)
    workB = check_geom_array_shape(geomB)
    if workB.shape != workA.shape:
        raise ValueError("Shape of input arrays does not match.")
    
    if prealign:
        workB = align(workA, workB)

    rmsd = np.linalg.norm(workA - workB)
    rmsd = rmsd / np.sqrt(workA.shape[0])
    return rmsd


def align(geomA, geomB):
    """ Align a geometry to match reference conformation.

    Args:
        geomA: Reference geometry.
        geomB: Geometry to align.

    Returns:
        geomB_rot: geomB aligned to minimize RMSD from geomA.
    """
    # Check array dimensions
    ndim = geomB.ndim
    workA = check_geom_array_shape(geomA)
    workB = check_geom_array_shape(geomB)
    if workB.shape != workA.shape:
        raise ValueError("Shape of input arrays does not match.")

    workA = recenter(workA)
    workB = recenter(workB)

    # Get rotation matrix
    R = kabsch(workA, workB)

    # Rotate second geometry
    geomB_rot = np.dot(workB, R)

    # Translation to geomA
    geomB_rot = geomB_rot + np.mean(geomA, axis=0)

    if ndim == 1:
        geomB_rot = geomB_rot.flatten()

    return geomB_rot


def recenter(geom):
    """ Recenter a set of points. """
    centroid = np.mean(geom, axis=0)
    return geom - centroid


def kabsch(geomA, geomB):
    """ Solve Wahba's problem using Kabsch algorithm.

    Args:
        geomA: First set of points, (npoint x dim) array.
        geomB: Second set of points, (npoint x dim) array.

    Returns:
        R: Optimal rotation matrix.
    """
    assert geomA.shape == geomB.shape

    # Translation
    Ac = recenter(geomA)
    Bc = recenter(geomB)

    # Covariance matrix
    cov_mat = np.dot(Ac.T, Bc)

    # Rotation
    U, S, Vt = np.linalg.svd(cov_mat)
    R = np.dot(Vt.T, U.T)

    # correct for possible reflection
    if np.linalg.det(R) < 0:
        Vt[-1,:] *= -1
        R = np.dot(Vt.T, U.T)

    return R 


def check_geom_array_shape(geom):
    """ Check dimensions of input array and reshape if needed.

    If a 1D array is passed reshape it to natom x 3 dimension. If 2D array is
    passed do nothing. Otherwise raise error.
    """
    if geom.ndim == 1 and geom.shape[0]%3 == 0:
        natom = geom.shape[0] // 3
        work = geom.reshape(natom, 3)
    elif geom.ndim == 2:
        work = geom
    else:
        raise ValueError("Bad shape of input array.")
    return work
