all:
	-latexmk test_results_document.tex -pdf
	-$(MAKE) -C slides
	mv test_results_document.pdf ../test_results_document.pdf 
	./copy_to_opendihu-ci.sh
	
continuous:
	latexmk test_results_document.tex -pvc -pdf -synctex=1 < /dev/null

clean:
	-rm -rf *.out *.log *.aux
	$(MAKE) clean -C slides
