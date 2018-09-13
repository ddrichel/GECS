gecs: gecs.cpp gecs.hpp rare.cpp alglib/alglibinternal.cpp alglib/alglibinternal.h alglib/ap.cpp alglib/ap.h alglib/specialfunctions.cpp alglib/specialfunctions.h alglib/stdafx.h
	g++ gecs.cpp -o gecs -O3 -lm -march=native

install: binge
.PHONY: install

