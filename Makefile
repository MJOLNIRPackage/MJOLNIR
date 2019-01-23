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


FILE := $(shell ls -t dist/* | head -1)
version: 
	make wheel
	git tag -a $(shell python cut.py $(shell ls -t dist/* | head -1))
	twine upload $(FILE) -r testpypi
	twine upload $(FILE) -r pypi


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


