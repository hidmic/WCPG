# coding: utf8

"""
Python wrapper for the Worst-Case Peak Gain library
"""

__author__ = "Thibault Hilaire, Anastasia Volkova"
__copyright__ = "Copyright 2015, FIPOgen Project, LIP6"
__credits__ = ["Thibault Hilaire", "Anastasia Volkova"]

__license__ = "GPL v3"
__version__ = "1.0"
__maintainer__ = "Thibault Hilaire"
__email__ = "thibault.hilaire@lip6.fr"
__status__ = "Beta"


import ctypes
import ctypes.util
from numpy import empty, float64, mat

# load the WCPG library, if installed
if not ctypes.util.find_library('wcpg'):
	raise ValueError("The WCPG library cannot be found (is it installed?")
_WCPGlib = ctypes.CDLL(ctypes.util.find_library('wcpg'))

# int WCPG_ABCD(double *W, double *A, double *B, double *C, double *D, uint64_t n, uint64_t p, uint64_t q)
_WCPGfunABCD = _WCPGlib.WCPG_ABCD
_WCPGfunABCD.argtypes = (5 * (ctypes.POINTER(ctypes.c_double),) + 3 * (ctypes.c_uint64,))
def WCPG_ABCD(A, B, C, D):
	"""Compute the WCPG from the matrices A, B, C, D
	A,B,C and D are numpy matrices or array of the right size"""
	# get the sizes
	n = A.shape[1]
	p,q = D.shape
	# get the pointer to the double arrays
	pA = A.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
	pB = B.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
	pC = C.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
	pD = D.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
	# run the function to fill the empty array W
	W = empty((p, q), dtype=float64)
	ret = _WCPGfunABCD(W.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), pA, pB, pC, pD, n, p, q)
	if ret == 0:
		raise ValueError("Something went wrong during the WCPG evaluation...")

	return mat(W)




# int WCPG_tf_initial(double *W, double *num, double *denum, uint64_t Nb, uint64_t Na)
_WCPGfunTF = _WCPGlib.WCPG_tf_initial
_WCPGfunTF.argtypes = (3 * (ctypes.POINTER(ctypes.c_double),) + 2 * (ctypes.c_uint64,))
def WCPG_TF(num, den):
	# TODO: finish and add the test function
	pass




