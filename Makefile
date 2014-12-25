include config.h

install: $(PREFIX)/lib
	cd src && $(MAKE) install

test::
	cd test && $(MAKE) test

$(PREFIX)/lib:
	/bin/mkdir $(PREFIX)/lib
