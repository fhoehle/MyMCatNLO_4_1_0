# Use this file as follows:
# gmake -f Makefile EXTRAOBJ=<alpha,linux>.o VPATH=<vpath> <EXENAME>
# See below for a list of name of executables. This is usually unnecessary,
# the relevant operations being done by the scripts. If done manually, the
# proper <vpath> must be entered at runtime
#
# Replace all instances of
#         g77 -w -fno-automatic
# with
#         gfortran -w -O2 -fno-automatic 
# in this file if g77 is not available or outdated

ifeq ($(shell uname),AIX)
F77=xlf -qextname -qflttrap=overflow:zerodivide:invalid:enable -O3 -qstrict \
#       -qautodbl=dblpad
SYSOBJ=
AUTODBL=-qautodbl=dblpad
endif
ifeq ($(shell uname),SunOS)
F77= f77 -fnonstd
SYSOBJ=
endif
ifeq ($(shell uname),Linux)
ifeq ($(COMPILERTYPE),xGFORTRAN)
# Use with 
F77= gfortran -w -O2 -fno-automatic 
# has been tested extensively with IPROC=-16XX and exact mass depedence.
#F77= g77 -w -fno-automatic
else
#F77= g77 -w -fno-automatic
F77= gfortran -w -O2 -fno-automatic
endif
CPP= g++ $(INCLOPTION)
CC= gcc
SYSOBJ=trapfpe.o
endif
ifeq ($(shell uname),HP-UX)
F77= g77 -w
SYSOBJ=
endif
ifeq ($(shell uname),OSF1)
F77= f77 
CPP= g++ $(INCLOPTION)
CC= gcc
SYSOBJ=
endif
ifeq ($(shell uname),Darwin)
ifeq ($(COMPILERTYPE),xGFORTRAN)
F77= gfortran -w -O2 -fno-automatic
else
F77= g77 -w -fno-automatic
endif
CPP= g++ $(INCLOPTION)
CC= gcc
endif

DEBUG=
FF=$(F77) $(DEBUG)

LIBS=`cernlib pdflib804 mathlib`
#$LIBS gets replaced at compilation time by
# /cern/pro/lib/libpdflib804.a /cern/pro/lib/libmathlib.a 
# /cern/pro/lib/libpacklib.a -L/usr/local/lib -lshift -lnsl -lcrypt -ldl
#on machines with a CERN installation. With a non-CERN installation,
#replace the definition above with 
#LIBS=/your/path/to/pdflib/libpdflib804.a .....

%.o: $(SRCDIR)/%.f
	$(F77) -I$(INCDIR) $(DEBUG) $(AUTODBL) -c $<
%.o: $(COMSRC)/%.f
	$(F77) -I$(INCDIR) $(DEBUG) $(AUTODBL) -c $<
%.o: $(ANADIR)/%.f
	$(F77) -I$(INCDIR) $(DEBUG) $(AUTODBL) -c $<
%.o: $(HWSDIR)/%.f
	$(F77) -I$(INCDIR) $(DEBUG) $(AUTODBL) -c $<

%.o: $(SRCDIR)/%.for
	$(F77) -I$(INCDIR) $(DEBUG) $(AUTODBL) -c $<
%.o: $(COMSRC)/%.for
	$(F77) -I$(INCDIR) $(DEBUG) $(AUTODBL) -c $<

%.o: $(SRCDIR)/%.cc
	$(CPP) -I$(INCDIR) $(DEBUG) -c $<
%.o: $(COMSRC)/%.cc
	$(CPP) -I$(INCDIR) $(DEBUG) -c $<
%.o: $(ANADIR)/%.cc
	$(CPP) -I$(INCDIR) $(DEBUG) -c $<

%.o: $(SRCDIR)/%.c
	$(CC) -I$(INCDIR) $(DEBUG) -c $^
%.o: $(COMSRC)/%.c
	$(CC) -I$(INCDIR) $(DEBUG) -c $^
%.o: $(ANADIR)/%.c
	$(CC) -I$(INCDIR) $(DEBUG) -c $<


VBFILES=mcatnlo_vbmain.o mcatnlo_vbxsec.o mcatnlo_vbpdks.o mcatnlo_helas2.o
QQFILES=mcatnlo_qqmain.o mcatnlo_qqxsec.o mcatnlo_helas2.o
HGFILES=mcatnlo_hgmain.o mcatnlo_hgxsec.o mcatnlo_chaplin_dummy.o
HGMFILES=mcatnlo_hgmain.o mcatnlo_hgxsec.o mcatnlo_chaplin11.o
SBFILES=mcatnlo_sbmain.o mcatnlo_sbxsec.o
LLFILES=mcatnlo_llmain.o mcatnlo_llxsec.o
VHFILES=mcatnlo_vhmain.o mcatnlo_vhxsec.o
STFILES=mcatnlo_stmain.o mcatnlo_stxsec.o mcatnlo_helas2.o
WTDRFILES=mcatnlo_wtmain_dr.o mcatnlo_wtxsec_dr.o mcatnlo_helas2.o
WTDSFILES=mcatnlo_wtmain_ds.o mcatnlo_wtxsec_ds.o mcatnlo_helas2.o
HTDRFILES=mcatnlo_htmain_dr.o mcatnlo_htxsec_dr.o mcatnlo_helas2.o
HTDSFILES=mcatnlo_htmain_ds.o mcatnlo_htxsec_ds.o mcatnlo_helas2.o
UTIFILES=mcatnlo_date.o mcatnlo_int.o mcatnlo_uxdate.o mcatnlo_uti.o \
         mcatnlo_str.o $(EXTRAOBJ)
LUTIFILES=mcatnlo_date.o mcatnlo_int.o mcatnlo_uxdate.o mcatnlo_uti.o \
         mcatnlo_str.o $(EXTRAOBJ)
PDFFILES=mcatnlo_pdftomlm.o mcatnlo_libofpdf.o dummies.o 
CPDFFILES=mcatnlo_mlmtopdf.o dummies.o 
LPDFFILES=mcatnlo_mlmtolha.o dummies.o 
HWFILES=mcatnlo_hwdriver.o mcatnlo_hwlhin.o \
        mcatnlo_str.o $(HWUTI)

QQNLO_EXE_THISLIB : $(QQFILES) $(UTIFILES) $(PDFFILES) $(SYSOBJ)
	$(FF) $^ $(EXTRAPATHS) $(EXTRALIBS) -o $@

QQNLO_EXE_PDFLIB : $(QQFILES) $(UTIFILES) $(CPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LIBS) $(EXTRAPATHS) $(EXTRALIBS) -o $@

QQNLO_EXE_LHAPDF : $(QQFILES) $(LUTIFILES) $(LPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LHALIB) $(EXTRAPATHS) $(EXTRALIBS) -o $@

VVNLO_EXE_THISLIB : $(VBFILES) $(UTIFILES) $(PDFFILES) $(SYSOBJ)
	$(FF) $^ $(EXTRAPATHS) $(EXTRALIBS) -o $@

VVNLO_EXE_PDFLIB : $(VBFILES) $(UTIFILES) $(CPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LIBS) $(EXTRAPATHS) $(EXTRALIBS) -o $@

VVNLO_EXE_LHAPDF : $(VBFILES) $(LUTIFILES) $(LPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LHALIB) $(EXTRAPATHS) $(EXTRALIBS) -o $@

HGNLO_EXE_THISLIB : $(HGFILES) $(UTIFILES) $(PDFFILES) $(SYSOBJ)
	$(FF) $^ $(EXTRAPATHS) $(EXTRALIBS) -o $@

HGNLO_EXE_PDFLIB : $(HGFILES) $(UTIFILES) $(CPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LIBS) $(EXTRAPATHS) $(EXTRALIBS) -o $@

HGNLO_EXE_LHAPDF : $(HGFILES) $(LUTIFILES) $(LPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LHALIB) $(EXTRAPATHS) $(EXTRALIBS) -o $@

HGMNLO_EXE_THISLIB : $(HGMFILES) $(UTIFILES) $(PDFFILES) $(SYSOBJ)
	$(FF) $^ $(EXTRAPATHS) $(EXTRALIBS) -o $@

HGMNLO_EXE_PDFLIB : $(HGMFILES) $(UTIFILES) $(CPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LIBS) $(EXTRAPATHS) $(EXTRALIBS) -o $@

HGMNLO_EXE_LHAPDF : $(HGMFILES) $(LUTIFILES) $(LPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LHALIB) $(EXTRAPATHS) $(EXTRALIBS) -o $@

SBNLO_EXE_THISLIB : $(SBFILES) $(UTIFILES) $(PDFFILES) $(SYSOBJ)
	$(FF) $^ $(EXTRAPATHS) $(EXTRALIBS) -o $@

SBNLO_EXE_PDFLIB : $(SBFILES) $(UTIFILES) $(CPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LIBS) $(EXTRAPATHS) $(EXTRALIBS) -o $@

SBNLO_EXE_LHAPDF : $(SBFILES) $(LUTIFILES) $(LPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LHALIB) $(EXTRAPATHS) $(EXTRALIBS) -o $@

LLNLO_EXE_THISLIB : $(LLFILES) $(UTIFILES) $(PDFFILES) $(SYSOBJ)
	$(FF) $^ $(EXTRAPATHS) $(EXTRALIBS) -o $@

LLNLO_EXE_PDFLIB : $(LLFILES) $(UTIFILES) $(CPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LIBS) $(EXTRAPATHS) $(EXTRALIBS) -o $@

LLNLO_EXE_LHAPDF : $(LLFILES) $(LUTIFILES) $(LPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LHALIB) $(EXTRAPATHS) $(EXTRALIBS) -o $@

VHNLO_EXE_THISLIB : $(VHFILES) $(UTIFILES) $(PDFFILES) $(SYSOBJ)
	$(FF) $^ $(EXTRAPATHS) $(EXTRALIBS) -o $@

VHNLO_EXE_PDFLIB : $(VHFILES) $(UTIFILES) $(CPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LIBS) $(EXTRAPATHS) $(EXTRALIBS) -o $@

VHNLO_EXE_LHAPDF : $(VHFILES) $(LUTIFILES) $(LPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LHALIB) $(EXTRAPATHS) $(EXTRALIBS) -o $@

STNLO_EXE_THISLIB : $(STFILES) $(UTIFILES) $(PDFFILES) $(SYSOBJ)
	$(FF) $^ $(EXTRAPATHS) $(EXTRALIBS) -o $@

STNLO_EXE_PDFLIB : $(STFILES) $(UTIFILES) $(CPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LIBS) $(EXTRAPATHS) $(EXTRALIBS) -o $@

STNLO_EXE_LHAPDF : $(STFILES) $(LUTIFILES) $(LPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LHALIB) $(EXTRAPATHS) $(EXTRALIBS) -o $@

WTDRNLO_EXE_THISLIB : $(WTDRFILES) $(UTIFILES) $(PDFFILES) $(SYSOBJ)
	$(FF) $^ $(EXTRAPATHS) $(EXTRALIBS) -o $@

WTDRNLO_EXE_PDFLIB : $(WTDRFILES) $(UTIFILES) $(CPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LIBS) $(EXTRAPATHS) $(EXTRALIBS) -o $@

WTDRNLO_EXE_LHAPDF : $(WTDRFILES) $(LUTIFILES) $(LPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LHALIB) $(EXTRAPATHS) $(EXTRALIBS) -o $@

WTDSNLO_EXE_THISLIB : $(WTDSFILES) $(UTIFILES) $(PDFFILES) $(SYSOBJ)
	$(FF) $^ $(EXTRAPATHS) $(EXTRALIBS) -o $@

WTDSNLO_EXE_PDFLIB : $(WTDSFILES) $(UTIFILES) $(CPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LIBS) $(EXTRAPATHS) $(EXTRALIBS) -o $@

WTDSNLO_EXE_LHAPDF : $(WTDSFILES) $(LUTIFILES) $(LPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LHALIB) $(EXTRAPATHS) $(EXTRALIBS) -o $@

HTDRNLO_EXE_THISLIB : $(HTDRFILES) $(UTIFILES) $(PDFFILES) $(SYSOBJ)
	$(FF) $^ $(EXTRAPATHS) $(EXTRALIBS) -o $@

HTDRNLO_EXE_PDFLIB : $(HTDRFILES) $(UTIFILES) $(CPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LIBS) $(EXTRAPATHS) $(EXTRALIBS) -o $@

HTDRNLO_EXE_LHAPDF : $(HTDRFILES) $(LUTIFILES) $(LPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LHALIB) $(EXTRAPATHS) $(EXTRALIBS) -o $@

HTDSNLO_EXE_THISLIB : $(HTDSFILES) $(UTIFILES) $(PDFFILES) $(SYSOBJ)
	$(FF) $^ $(EXTRAPATHS) $(EXTRALIBS) -o $@

HTDSNLO_EXE_PDFLIB : $(HTDSFILES) $(UTIFILES) $(CPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LIBS) $(EXTRAPATHS) $(EXTRALIBS) -o $@

HTDSNLO_EXE_LHAPDF : $(HTDSFILES) $(LUTIFILES) $(LPDFFILES) $(SYSOBJ)
	$(FF) $^ $(LHALIB) $(EXTRAPATHS) $(EXTRALIBS) -o $@

MC_EXE_THISLIB : $(HWFILES) $(HERWIGVER) $(PDFFILES)
	$(FF) $^ $(EXTRAPATHS) $(EXTRALIBS) -o $@

MC_EXE_PDFLIB : $(HWFILES) $(HERWIGVER) $(CPDFFILES) 
	$(FF) $^ $(LIBS) $(EXTRAPATHS) $(EXTRALIBS) -o $@

MC_EXE_LHAPDF : $(HWFILES) $(HERWIGVER) $(LPDFFILES)  
	$(FF) $^ $(LHALIB) $(EXTRAPATHS) $(EXTRALIBS) -o $@