plaquette:
	python iqcd.py code +cold +loop n=100 { +heatbath +plaquette }
	python iqcd.py compile
	./latest.exe > latest.log	
	python ibootstrap.py latest.log '"plaquette"'
	python iplot.py -p quartz

topological:
	python iqcd.py code +hot +heatbath steps=100 +loop n=100 { +ape_smear steps=1 +topological_charge }
	python iqcd.py compile
	./latest.exe > latest.log
	python ibootstrap.py latest.log '"topological_charge"'
	python isee.py topological_charge_00000.vtk

ising3d:
	./ising3d.exe > ising3d.log &
	python isee.py ising3d.vtk

