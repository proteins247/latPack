

################################################################
# doxygen support
################################################################

include $(top_srcdir)/aminclude.am

MOSTLYCLEANFILES = $(DX_CLEANFILES)

################################################################
# new doc rule to ensure that the picture is copied
################################################################

doc: doxygen-doc
if DX_COND_html
	cp $(top_srcdir)/misc/doxygen/bioinf-fr-logo.png $(top_srcdir)/$(DX_DOCDIR)/html/.
endif

################################################################
# pkg-config
################################################################

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = biu.pc


DOXYGENSUPPORTFILES = acinclude.m4 aminclude.am

EXTRA_DIST = $(DOXYGENSUPPORTFILES)

tests:	all
	@echo "-------------------------------------------------------------------"
	@echo "Run test"
	@make tests -s -C tests
	@echo "-------------------------------------------------------------------"

test: tests

clean-local:
	$(RM) $(top_srcdir)/src/*.gc* $(top_srcdir)/tests/*.gc* $(MOSTLYCLEANFILES)

SUBDIRS = src misc tests
	