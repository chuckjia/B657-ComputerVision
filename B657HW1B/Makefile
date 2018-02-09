all: detect

detect: detect.cpp
	g++ detect.cpp -o detect -lpng -I.

clean:
	rm detect

test: detect.cpp
	g++ detect.cpp -o detect -lpng -I.
	./detect ic_6.png
	python2.7 a1eval.py ic_6.txt detected.txt
	
