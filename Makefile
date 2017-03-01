SHELL = /bin/bash -e 

all: build install

build:
	python setup.py build --executable="/usr/bin/env python"

bdist:
	python setup.py build --executable="/usr/bin/env python"
	python setup.py bdist --formats=egg

install:
	python setup.py install

develop:
	python setup.py develop

test: unittests cramtests

unittests:
	# Unit tests
	find tests/unit -name "*.py" | xargs nosetests

cramtests:
	# End-to-end tests
	find tests/cram -name "*.t" | xargs cram
	find tests/cram/test_ice_entries -name "*.t" | xargs cram

doc:
	sphinx-apidoc -T -f -o doc src/ && cd doc && make html
docs: doc

clean: doc-clean
	rm -rf dist/ build/ *.egg-info
	rm -rf doc/_build
	find . -name "*.pyc" | xargs rm -f
	rm -rf dist/
	rm -f nostests.xml
	rm -f pbtranscript/collapsing/C/c_branch.cpp
	rm -f pbtranscript/collapsing/C/intersection.cpp
	rm -f pbtranscript/collapsing/C/intersection_unique.cpp
	rm -f pbtranscript/io/C/SAMReaders.cpp

doc-clean:
	rm -f doc/*.html

pip-install:
	@which pip > /dev/null
	@pip freeze|grep 'pbtranscript=='>/dev/null \
      && ( pip uninstall -y pbtranscript \
        || pip uninstall -y pbtools.pbtranscript ) \
      || true
	@python setup.py build
	@pip install --no-index \
          --install-option="--install-scripts=$(PREFIX)/bin" \
          ./

.PHONY: all build bdist install develop test doc clean pip-install
