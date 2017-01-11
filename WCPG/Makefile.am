#######################################
# The list of libraries we are building seperated by spaces.
# The 'lib_' indicates that these build products will be installed
# in the $(libdir) directory. For example /usr/lib

AM_LDFLAGS = -lm -lgmp -lmpfr -lmpfi $(LAPACK_LIBS) $(BLAS_LIBS)

lib_LTLIBRARIES = libwcpg.la

#######################################
# Build information for each library

# Sources for WCPG
libwcpg_la_SOURCES = wcpg.c \
			   _wcpg.h \
			   aux_funcs.c \
			   clapack_linalg.c \
			   clapack_linalg.h \
			   glue_conversions.c \
			   glue_conversions.h \
			   inclusion_verif.c \
			   inclusion_verif.h \
			   mpfi_matrix.c \
			   mpfi_matrix.h \
			   mpfr_linalg.c \
			   mpfr_linalg.h \
			   mpfr_matrix.c \
			   clapack_functions_config.c  \
			   clapack_functions_config.h \
			   mpfr_eigendecomposition.c \
			   mpfr_eigendecomposition.h

# Linker options libTestProgram
libwcpg_la_LDFLAGS = -version-info 0:9:0 $(AM_LDFLAGS)

# Compiler options. Here we are adding the include directory
# to be searched for headers included in the source code.
libwcpg_la_CPPFLAGS = -I$(top_srcdir)/include

LIBTOOL_DEPS = @LIBTOOL_DEPS@
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status libtool