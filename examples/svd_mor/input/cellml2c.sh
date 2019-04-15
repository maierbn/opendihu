file=$1
/store/software/opencor/OpenCOR-0-5-Linux/bin/OpenCOR -c CellMLTools::export $file /store/software/opencor/OpenCOR-0-5-Linux/formats/C.xml > ${file%.*}.c
