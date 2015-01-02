include config.h

install: $(PREFIX)/lib
	cd src && $(MAKE) install

clean:
	-cd src && $(MAKE) clean
	-cd test && $(MAKE) clean

test::
	cd test && $(MAKE) test

docs::
	doxygen doxyfile.in

$(PREFIX)/lib:
	/bin/mkdir $(PREFIX)/lib
