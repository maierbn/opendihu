all:
	latexmk doc.tex -pdf
continuous:
	latexmk doc.tex -pdf -pvc < /dev/null
clean:
	rm -f doc.aux doc.bbl doc.blg doc.fdb_latexmk doc.fls doc.log doc.out doc.synctex.gz 'doc.synctex.gz(busy)' doc.toc
