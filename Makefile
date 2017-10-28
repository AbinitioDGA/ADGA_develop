all: binR abinitiodga vq make_config 

make_config:
	cd make_configs/;./make_config_auto.sh

binR:
	if [ ! -d bin ] ; then mkdir bin ; fi

vq: make_config binR
	cd vq/; make

abinitiodga: make_config
	cd src/; make

test: all
	cd testsuit/; make clean; make

clean: make_config
	cd testsuit/; make clean
	cd src/; make clean
	cd vq/; make clean

pristine: make_config
	cd testsuit/; make pristine
	cd src/; make pristine
	cd vq/; make pristine
	rmdir bin

