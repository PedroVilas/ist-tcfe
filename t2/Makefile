#!/bin/bash
# ------------------------------------------------------------
# type "make" command in the Linux terminal create report.pdf
#
# type "make clean" to delete all generated files
# ------------------------------------------------------------

all:
	make -C mat
	make -C sim
	mv sim/op_tab.tex doc/op_tab.tex
	mv sim/opeq_tab.tex doc/opeq_tab.tex
	make -C doc
	cp doc/report.pdf .

clean:
	@rm -f data.txt
	make -C mat clean
	make -C sim clean
	make -C doc clean

.PHONY: all clean

