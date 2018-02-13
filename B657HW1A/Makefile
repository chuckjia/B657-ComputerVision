all: watermark 

watermark: watermark.cpp
	g++ watermark.cpp -o watermark -O3 -lpng -I.

clean:
	rm watermark 

addcheck: watermark.cpp
	g++ watermark.cpp -o watermark -O3 -lpng -I.
	number=1 ; while [[ $$number -le 10000 ]] ; do \
        	echo $$number ; \
		./watermark 1.3 sea.png sea_marked.png add $$number ; \
		./watermark 1.3 sea_marked.png no_file.png check $$number ; \
        	((number = number + 100)) ; \
    	done

checkother: watermark.cpp
	g++ watermark.cpp -o watermark -O3 -lpng -I.
	./watermark 1.3 sea.png sea_marked.png add 2; \
	number=1 ; while [[ $$number -le 10000 ]] ; do \
        	echo Seed = $$number ; \
		./watermark 1.3 sea_marked.png no_file.png check $$number ; \
        	((number = number + 100)) ; \
    	done
