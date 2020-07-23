PYTHON=python
path=examples
recursive=True

make:
	@echo Installing pygibbs...
	${PYTHON} setup.py install

clean-local:
	@echo removing local compiled files
	rm pygibbs/*.c pygl/*.html pygl/*.cpp

clean:
	@echo removing all compiled files
	${PYTHON} setup.py clean
	rm pygibbs/*.c pygl/*.html 
	
env:
	@echo creating conda environment...
	conda env create --file environment.yml
	# conda activate pygibbs
	@echo use make to install pygibbs

test:
	@echo testing pygibbs...
	cd pygibbs/tests/ && python installTests.py

nbtest:
	@echo testing example notebooks...
	@echo test $(path)
	cd pygibbs/tests/ && python testNotebooks.py --path $(path) --recursive $(recursive)


pypitest:
	@echo testing pystokes...
	python setup.py sdist bdist_wheel
	python -m twine upload --repository testpypi dist/*

pypi:
	@echo testing pystokes...
	python setup.py sdist bdist_wheel	
	python -m twine upload dist/*

