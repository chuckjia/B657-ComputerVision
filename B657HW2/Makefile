all: CImg.h a2.cpp
	g++ a2.cpp -Wextra -o a2 -I/opt/X11/include -L/opt/X11/lib -lX11 -lpthread -I. -Isiftpp -O3 siftpp/sift.cpp -ljpeg
	./a2

clean:
	rm a2