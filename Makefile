all: binR abinitiodga make_config 

make_config:
	cd make_configs/;./make_config_auto.sh

binR:
	if [ ! -d bin ] ; then mkdir bin ; fi

abinitiodga: make_config
	cd src/; make

test: all
	cd testsuit/; make clean; make

clean: make_config
	cd testsuit/; make clean
	cd src/; make clean

pristine: make_config
	cd testsuit/; make pristine
	cd src/; make pristine
	rmdir bin
