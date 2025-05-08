PYTHON=python
path=examples
recursive=True

make:
	@echo Installing pygl...
	${PYTHON} setup.py install

clean-local:
	@echo removing local compiled files
	rm pygl/*.c pygl/*.html pygl/*.cpp

clean:
	@echo removing all compiled files
	${PYTHON} setup.py clean
	rm pygl/*.c pygl/*.html 
	
env:
	@echo creating conda environment...
	conda env create --file environment.yml
	# conda activate pygl
	@echo use make to install pygl

test:
	@echo testing pygl...
	cd pygl/tests/ && python installTests.py

nbtest:
	@echo testing example notebooks...
	@echo test $(path)
	cd pygl/tests/ && python testNotebooks.py


pypitest:
	@echo testing pystokes...
	python setup.py sdist bdist_wheel
	python -m twine upload --repository testpypi dist/*

pypi:
	@echo testing pystokes...
	python setup.py sdist bdist_wheel	
	python -m twine upload dist/*

