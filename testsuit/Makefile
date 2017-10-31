test: tests.yaml runcommand
	./beetest

extendedtest: tests_extended.yaml runcommand
	./beetest tests_extended.yaml

alltests: tests.yaml tests_extended.yaml runcommand
	./beetest tests.yaml tests_extended.yaml 

runcommand:
	cp runcommand.default runcommand

clean: 
	make pristine

pristine: 
	cd Tests/oneband_eom; make clean
	cd Tests/oneband_eom_large; make clean
	cd Tests/oneband_chi; make clean
	cd Tests/oneband_eomchi; make clean
	cd Tests/oneband_vq_eomchi; make clean
	cd Tests/oneband_vq_chipath; make clean
	cd Tests/SVO_eom; make clean
	cd Tests/SVO_chi; make clean
	cd Tests/SVO_eomchi; make clean
	cd Tests/SVO_vq_eomchi; make clean

