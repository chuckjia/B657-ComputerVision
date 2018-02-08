all: detect

detect: detect.cpp
	g++ detect.cpp -o detect -lpng -I.

clean:
	rm detect

test: detect.cpp
	g++ detect.cpp -o detect -lpng -I.
	./detect ic_1.png
	./detect ic_2.png
	./detect ic_3.png
	./detect ic_4.png
	./detect ic_5.png
	./detect ic_6.png
	./detect ic_7.png
	./detect ic_8.png
	./detect ic_9.png

test1: detect.cpp
	g++ detect.cpp -o detect -lpng -I.
	./detect ic_1.png
	
