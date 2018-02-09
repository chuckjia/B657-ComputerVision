all: watermark 

watermark: watermark.cpp
	g++ watermark.cpp -o watermark -O3 -lpng -I.

clean:
	rm watermark 

testadd: watermark.cpp
	g++ watermark.cpp -o watermark -O3 -lpng -I.
	number=1 ; while [[ $$number -le 10000 ]] ; do \
        	echo $$number ; \
		./watermark 1.3 iu.png iu_marked.png add $$number ; \
		./watermark 1.3 iu_marked.png iu_marked.png check $$number ; \
        	((number = number + 100)) ; \
    	done

testcheckother: watermark.cpp
	g++ watermark.cpp -o watermark -O3 -lpng -I.
	./watermark 1.3 iu.png iu_marked.png add 2; \
	number=1 ; while [[ $$number -le 10000 ]] ; do \
        	echo Seed = $$number ; \
		./watermark 1.3 iu_marked.png no_file.png check $$number ; \
        	((number = number + 100)) ; \
    	done
	

testnoadd: watermark.cpp
	g++ watermark.cpp -o watermark -O3 -lpng -I.
	./watermark 1.3 iu.png no_file.png check 90

testoneadd: watermark.cpp
	g++ watermark.cpp -o watermark -O3 -lpng -I.
	./watermark 1.3 desert.png desert_marked.png add 90
	./watermark 1.3 desert.png desert_marked.png check 90
