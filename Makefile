.PHONY: help all install test clean doc uninstall wheel deploy deploy_test

all: install test

help:
	@echo "make all"
	@echo "    install Biofrost nd dependencies into current python environment"
	@echo "make install"
	@echo "    install Biofrost"
	@echo "make test"
	@echo "    unit test, run after installation"
	@echo "make clean"
	@echo "    clean python cache files"
	@echo "make doc"
	@echo "    build API documentation"

install:
	python3 setup.py install

test:
	python3 setup.py test

clean:
	python3 setup.py clean --all

doc:
	pip install -r docs/source/requirements.txt
	cd docs && sphinx-apidoc -e -M -f -o ./source ../biofrost
	cd docs && make html

uninstall:
	pip uninstall biofrost

wheel:
	make clean
	pip install wheel setuptools
	pip install twine
	python setup.py sdist bdist_wheel

deploy_test:
	twine check dist/*
	twine upload --repository-url https://test.pypi.org/legacy/ dist/*

deploy:
	twine check dist/*
	twine upload dist/*
