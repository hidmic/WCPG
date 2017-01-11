#######################################
# The list of executables we are building seperated by spaces
# the 'bin_' indicates that these build products will be installed
# in the $(bindir) directory. For test /usr/bin
#bin_PROGRAMS=testProgram

# Because a.out is only a sample program we don't want it to be installed.
# The 'noinst_' prefix indicates that the following targets are not to be
# installed.
noinst_PROGRAMS=testProgram

#######################################
# Build information for each executable. The variable name is derived
# by use the name of the executable with each non alpha-numeric character is
# replaced by '_'. So a.out becomes a_out and the appropriate suffex added.
# '_SOURCES' for test.

ACLOCAL_AMFLAGS=-I ../m4

# Sources for the a.out 
testProgram_SOURCES= testProgram.c 

# Libraries for a.out
testProgram_LDADD = $(top_srcdir)/WCPG/libwcpg.la

# Linker options for a.out
testProgram_LDFLAGS = -rpath `cd $(top_srcdir);pwd`/WCPG/.libs

# Compiler options for a.out
testProgram_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/WCPG -g -O0
