include config.h

install: lib
	cd src && $(MAKE) install

test::
	cd test && $(MAKE) test

lib:
	/bin/mkdir lib
