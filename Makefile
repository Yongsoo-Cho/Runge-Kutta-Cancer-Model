project: helix_break.o radiation.o tumour_growth.o radiation.h tumour_growth.h
	c++ helix_break.o radiation.o tumour_growth.o -o project -ltrapfpe -lpgplot -lcpgplot -lX11 -lm 

helix_break.o: helix_break.cpp
	c++ helix_break.cpp -c

radiation.o: radiation.cpp
	c++ radiation.cpp -c

tumour_growth.o: tumour_growth.cpp
	c++ tumour_growth.cpp -c





