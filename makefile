PYTHON=python
path=examples
recursive=True

make:
	@echo Installing pymaft...
	${PYTHON} setup.py install

clean-local:
	@echo removing local compiled files
	rm pymaft/*.c pymaft/*.html pymaft/*.cpp

clean:
	@echo removing all compiled files
	${PYTHON} setup.py clean
	rm pymaft/*.c pymaft/*.html 
	
env:
	@echo creating conda environment...
	conda env create --file environment.yml
	# conda activate pymaft
	@echo use make to install pymaft

test:
	@echo testing pymaft...
	cd pymaft && python installTests.py

nbtest:
	@echo testing example notebooks...
	@echo test $(path)
	cd examples && python testNotebooks.py --path $(path) --recursive $(recursive)
