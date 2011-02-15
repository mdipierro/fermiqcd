all:
	cd Examples; make all
	cd Converters; make all
backup:
	cd ../..; tar cvf fermiqcd_3.3.1.tar FermiQCD/Version_3.3.1/*; gzip fermiqcd_3.3.1.tar
clean:
	rm */*.exe */*.o */*~
launchpad:
	bzr push bzr+ssh://mdipierro@bazaar.launchpad.net/~mdipierro/qcd/development --use-existing-dir