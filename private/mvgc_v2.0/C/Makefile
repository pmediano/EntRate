# Makefile for MVGC mex files. Run 'make' from the command line
# in this directory. If build is successful, run 'mextest' from
# the Matlab prompt.
#
# Tested on gcc (Linux 64-bit) and mingw-w64 (Windows 64-bit);
# note that MSVC (still) does not fully support the C99
# standard, so you may have problems building with it. Should
# also build okay on MacOS with gcc (not sure about clang/LLVM).

ifdef ComSpec
	OS = WIN
else
	UNAME = $(shell uname)
	ifeq ($(UNAME),Linux)
		OS = LINUX
	else
		ifeq ($(UNAME),Darwin)
			OS = MAC
		else
			OS = UNKNOWN
		endif
	endif
endif

MEXEXT =  $(shell mexext)
MEXDIR = ../mex

MSRC   = findin_mex.c mvfilter_mex.c
MBIN   = $(patsubst %.c,$(MEXDIR)/%.$(MEXEXT),$(MSRC))

MFLAGS = -std=c99 -march=native
WFLAGS = -Wall -Wextra -Wconversion -Winline -pedantic-errors
DFLAGS = -D_POSIX_C_SOURCE=199309L -D_BSD_SOURCE -D$(OS)
ifeq ($(OS),WIN)
	OFLAGS    = -O3
	MEXCFLAGS = COMPFLAGS
	RM        = del /F /Q
	RMBIN     = $(subst /,\,$(MBIN))
else
	OFLAGS    = -O3 -flto
	MEXCFLAGS = CFLAGS
	RM        = rm -f
	RMBIN     = $(MBIN)
endif
GCCFLAGS = $(MFLAGS) $(OFLAGS) $(WFLAGS) $(DFLAGS)
MEXFLAGS = -O -largeArrayDims

CC  = gcc
MEX = mex

.PHONY: all clean diag

all: $(MBIN)

clean:
	$(RM) $(RMBIN)

diag:
	@echo "*** MSRC = " $(MSRC)
	@echo "*** MBIN = " $(MBIN)

$(MEXDIR)/%.$(MEXEXT): %.c
	@echo
	@echo Building $@
	$(MEX) $(MEXFLAGS) $(MEXCFLAGS)="\$$CFLAGS $(GCCFLAGS) -static" -outdir $(MEXDIR) $<
	@echo If build succeeded, you should run \'mextest\' in Matlab
