all: watermark 

watermark: watermark.cpp
	g++ watermark.cpp -o watermark -O3 -lpng -I.

clean:
	rm watermark 

testadd: watermark.cpp
	g++ watermark.cpp -o watermark -O3 -lpng -I.
	number = 1 ; while [[ $$number -le 100 ]] ; do \
        	echo $$number ; \
		./watermark 1.3 sea.png out.png add $$number ; \
		./watermark 1.1 out.png out_mag.png ; \
		./watermark 1.3 out.png output.png check $$number ; \
        	((number = number + 1)) ; \
    	done
	

testnoadd: watermark.cpp
	g++ watermark.cpp -o watermark -O3 -lpng -I.
	./watermark 1.3 sea.png output.png check 90
