CLEANFILES = out out.[1-9]*  

test:
	make clean; cd ../../../; make abinitiodga
	./run.sh

clean:
	rm -f $(CLEANFILES)

