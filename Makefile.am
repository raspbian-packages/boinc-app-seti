#  $Id: Makefile.am,v 1.2.2.1 2005/10/26 18:23:03 korpela Exp $
#
#

AUTOMAKE_OPTIONS = foreign

## make sure rebuilding uses the right aclocal-switches
ACLOCAL_AMFLAGS = -I m4

if ENABLE_CLIENT
  CLIENT_SUBDIRS = client
if ENABLE_TESTS
  CLIENT_SUBDIRS +=tools
endif
endif

if ENABLE_SERVER
  SERVER_SUBDIRS = splitter_pfb assimilator validate
endif

if USE_MY_IMAGELIBS
  IMDIR = image_libs
endif

if USE_MY_GLUT
  GLUTDIR = glut
endif


SUBDIRS = $(IMDIR) $(CLIENT_SUBDIRS) $(SERVER_SUBDIRS)

dist-hook:
	rm -rf `find $(distdir) -name .backup -prune -o -name '*~' -o -name CVS`

client-bin:
	cd client && make client-bin client-bin-gz

EXTRA_DIST =                                    \
	glut					\
	image_libs				\
	db					\
	tools					\
        client/win_build                                 \
        client/test_workunits                                   \
    INSTALL

