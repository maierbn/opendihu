scp ../test_results_document.pdf opendihu-ci:software/apache/public-html/documents
scp ../test_results_slides.pdf opendihu-ci:software/apache/public-html/documents
echo "Document and slides created on $(date '+%d/%m/%Y %H:%M:%S').<br>" > date.txt
scp date.txt opendihu-ci:software/apache/public-html/documents

