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
	cd Tests/oneband; make clean
	cd Tests/onebandVq; make clean
	cd Tests/SVO; make clean

