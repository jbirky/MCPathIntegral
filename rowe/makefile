CC=g++

# LFLAGS = -lm /usr/share/doc/blas -lcblas

LFLAGS = -lm -L/System/Library/Frameworks/Python.framework/Versions/2.6/Extras/lib/python/scipy/lib/blas/ -lcblas 

QMC_Psi: QMC_Psi.cpp
	$(CC) -o QMC_Psi QMC_Psi.cpp $(LFLAGS); ./QMC_Psi
	