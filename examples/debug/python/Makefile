
basedir := /store/software/opendihu
all:
	rm test && g++ -isystem $(basedir)/dependencies/python/install/include/python3.6m test.cpp -L$(basedir)/dependencies/python/install/lib -lpython3.6m -lpthread -ldl -lutil -o test && ./test
