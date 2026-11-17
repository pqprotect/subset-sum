PYTHON ?= python3
PDFLATEX ?= pdflatex
BIBTEX ?= bibtex

.PHONY: paper verify clean

paper:
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error main.tex
	$(BIBTEX) main
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error main.tex
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error main.tex
	cp main.pdf preprint.pdf

verify:
	$(PYTHON) verify_classical.py
	$(PYTHON) verify_qracm.py
	$(PYTHON) verify_qraqm.py

clean:
	rm -f main.aux main.blg main.log main.out main.pdf
