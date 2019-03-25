# Minimal makefile for Sphinx documentation
#

test:
	python -m pytest -vv MJOLNIR

quicktest:
	python -m pytest -vv MJOLNIR --quick

tutorials:
	python clean.py docs/Tutorials/Advanced
	python clean.py docs/Tutorials/Instrument
	python clean.py docs/Tutorials/Quick
	python clean.py docs/Tutorials/Tools
	./Tutorials/tutorials
	make html

fulltest:
	python -m pytest -vv MJOLNIR


wheel:
	python setup.py sdist

version = $(shell python cut.py $(shell ls -t dist/* | head -1))

upload:
	twine upload $(shell ls -t dist/* | head -1) -r testpypi
	twine upload $(shell ls -t dist/* | head -1) -r pypi


version: 
	echo 'Creating version $(version)'
	python Update.py $(version)
	make tutorials
	git add setup.py docs/conf.py docs/Tutorials/* docs/index.rst
	git commit -m 'Update version'
	git tag -a $(version) -m \'$(version)\'
	make wheel
	git push
	git push --tags


# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = MJOLNIR
SOURCEDIR     = docs
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)


