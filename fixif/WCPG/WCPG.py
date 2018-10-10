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


class wcpg_result(ctypes.Structure):
	_fields_ = [('N', ctypes.c_int),                    # number of iterations
		('one_minus_rhoA', ctypes.c_double),    # 1 - rho(A) where rho(A) is the spectral radius of A
		('maxSN', ctypes.c_double),             # max(S_N)
	    ('minSN', ctypes.c_double),             # min(S_N)
		('time_overall', ctypes.c_double),       # Overall time spent
		('time_Ncomp', ctypes.c_double),         # computation time
		('time_Summation', ctypes.c_double),     # Summation time
		('inversion_Iter', ctypes.c_int),     # Inversion iterations
		('maxprec_PN', ctypes.c_uint64),         # Maximum precision of U = inv(V)
		('maxprec_U', ctypes.c_uint64),         # Maximum precision of P_N
		('maxprec_SN', ctypes.c_uint64)]        # Maximum precision of S_N

	def getdict(self):
		return dict((field, getattr(self, field)) for field, _ in self._fields_)


# load the WCPG library, if installed
if not ctypes.util.find_library('wcpg'):
	raise ValueError("The WCPG library cannot be found (is it installed?")
_WCPGlib = ctypes.CDLL(ctypes.util.find_library('wcpg'))


# -- WCPG_ABCD --
# compute the WCPG from state-space matrices
# int WCPG_ABCD(double *W, double *A, double *B, double *C, double *D, uint64_t n, uint64_t p, uint64_t q);
# int WCPG_ABCD_res(double *W, double *A, double *B, double *C, double *D, uint64_t n, uint64_t p, uint64_t q, wcpg_result_out* res);
_WCPGfunABCD = _WCPGlib.WCPG_ABCD
_WCPGfunABCD.argtypes = (5 * (ctypes.POINTER(ctypes.c_double),) + 3 * (ctypes.c_uint64,))
_WCPGfunABCD_res = _WCPGlib.WCPG_ABCD_res
_WCPGfunABCD_res.argtypes = (5 * (ctypes.POINTER(ctypes.c_double),) + 3 * (ctypes.c_uint64,) + (ctypes.POINTER(wcpg_result),))


def WCPG_ABCD(A, B, C, D, result_info=None):
	"""Compute the WCPG from the matrices A, B, C, D
	A,B,C and D are numpy matrices or array of the right size
	if res is given, this dictionary is filled with result information"""
	# get the sizes
	n = A.shape[1]
	p, q = D.shape
	# get the pointer to the double arrays
	pA = A.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
	pB = B.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
	pC = C.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
	pD = D.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
	# get the result
	res = wcpg_result()
	# run the function to fill the empty array W
	W = empty((p, q), dtype=float64)
	if result_info is None:
		ret = _WCPGfunABCD(W.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), pA, pB, pC, pD, n, p, q)
	else:
		ret = _WCPGfunABCD_res(W.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), pA, pB, pC, pD, n, p, q, ctypes.pointer(res))
		result_info.update(res.getdict())
	if ret == 0:
		raise ValueError("Something went wrong during the WCPG evaluation...")

	return mat(W)



# -- WCPG_TF --
# compute the WCPG from transfer function coefficients
# int WCPG_tf_initial(double *W, double *num, double *denum, uint64_t Nb, uint64_t Na)
_WCPGfunTF = _WCPGlib.WCPG_tf
_WCPGfunTF.argtypes = (3 * (ctypes.POINTER(ctypes.c_double),) + 3 * (ctypes.c_uint64,))


def WCPG_TF(num, den):
	"""Compute the WCPG from the numerator and denominator of the transfer function
	num and den are numpy vectors or array of the right size"""
	# get the sizes
	Na = den.size
	Nb = num.size

	# get the pointer to the double arrays
	pNum = num.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
	pDen = den.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
	# run the function to fill the empty array W
	W = empty((1, 1), dtype=float64)
	ret = _WCPGfunTF(W.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), pNum, pDen, Nb, Na, 1)
	if ret == 0:
		raise ValueError("Something went wrong during the WCPG evaluation...")

	return mat(W)




