all: CImg.h a2.cpp
	g++ a2.cpp -Wextra -o a2 -I/opt/X11/include -L/opt/X11/lib -lX11 -lpthread -I. -Isiftpp -O3 siftpp/sift.cpp -ljpeg
	./a2

fast: CImg.h a2.cpp
	g++ a2.cpp -Wextra -o a2 -I/opt/X11/include -L/opt/X11/lib -lX11 -lpthread -I. -Isiftpp -O3 siftpp/sift.cpp -ljpeg -include CImg.h -static 

clean:
	rm a2

header:
	g++ CImg.h -I/opt/X11/include -L/opt/X11/lib 