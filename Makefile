.PHONY: help all install test clean doc uninstall

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
	cd docs && sphinx-apidoc -f -o ./source ../biofrost
	cd docs && make html

uninstall:
	pip uninstall biofrost