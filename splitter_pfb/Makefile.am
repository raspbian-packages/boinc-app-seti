CC=gcc
BOINCDIR=@BOINCDIR@
INFORMIXDIR=@INFORMIXDIR@
SETILIB_PATH=@SETILIB_PATH@
SETILIB_LIBS=@SETILIB_LIBS@
SETIHOME=..

LINKOPTIONS=

GSL_LIBS = -lgsl -lgslcblas -lgsl

DBLIBS=@INFORMIX_LIBS@ -lm -lstdc++


LINKOPTIONS=-Xlinker -R -Xlinker $(INFORMIXDIR)/lib:$(INFORMIXDIR)/lib/esql:$(LD_LIBRARY_PATH)


AM_CFLAGS= -g -O3 -Wall $(INCLUDE_DIRS) -DUSE_INFORMIX @PTHREAD_CFLAGS@ -ISETILIB_PATH/include


SYSLIBS = -lcrypto -ldl

BOINCLIBS= -L$(BOINCDIR)/sched -lsched @MYSQL_LIBS@ @PTHREAD_LIBS@ -L$(BOINCDIR)/lib -lboinc_crypt -lboinc -L$(SSLDIR) -lcrypto -lssl

SPLITLIBS= $(SETILIB_LIBS) \
        $(DBLIBS) \
	$(BOINCLIBS) \
        $(SYSLIBS) \
	-lchealpix \
	-lfftw3f

noinst_PROGRAMS = mb_splitter


mb_splitter_SOURCES=mb_angdist.cpp \
		    mb_message.cpp \
		    mb_splitter.cpp \
		    mb_read_blocks_dr2.cpp \
		    mb_read_blocks_guppi.cpp \
		    mb_wufiles.cpp \
		    mb_dotransform.cpp \
		    mb_validrun.cpp \
		    mb_coords.cpp \
		    cmap_interp.cpp \
		    ../db/schema_master.cpp \
	            ../db/sqlifx.cpp \
		    ../db/sqlrow.cpp \
	            ../db/sqlblob.cpp ../db/sqlint8.cpp \
	            ../db/xml_util.cpp \
	            ../db/app_config.cpp  \
	             ../client/seti_header.cpp \
	             ../client/timecvt.cpp \
	             ../client/lcgamm.cpp \
		     ../client/hr_min_sec.o

mb_splitter_CXXFLAGS= \
	   	-O3 \
       	-I$(SETIHOME) -I$(SETIHOME)/client \
		-I$(SETIHOME)/db @MYSQL_CFLAGS@ \
		-I$(HEALPIX)/include @INFORMIX_CFLAGS@ \
		@SETILIB_CFLAGS@ @BOINC_CFLAGS@ \
		-I$(BOINCDIR)/tools -I$(BOINCDIR)/sched \
		-I$(BOINCDIR)/db -I$(CFITSIO) \
		-I$(SETI_GBT)/src
mb_splitter_CFLAGS=$(mb_splitter_CXXFLAGS)
mb_splitter_LDFLAGS=-static $(AM_LDFLAGS) -L$(HEALPIX)/lib -L$(SETI_GBT)/lib -lsetigbt -lsla -L$(CFITSIO) -lcfitsio $(LINKOPTIONS)
mb_splitter_LDADD=$(SPLITLIBS) $(GSL_LIBS)

../db/sqlifx.cpp: ../db/sqlifx.ec
	$(INFORMIXDIR)/bin/esql -e $<
	mv sqlifx.c $*.cpp

