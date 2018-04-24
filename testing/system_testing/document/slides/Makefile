all:
	latexmk presentation.tex -pdf
	mv presentation.pdf ../../test_results_slides.pdf 
clean:
	rm *.aux *.fdb_latexmk *.log *.out
continuous:
	latexmk presentation.tex -pvc -pdf -synctex=1 < /dev/null
