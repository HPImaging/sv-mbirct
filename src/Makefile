CC = icc

#CFLAGS := -Wall -qopenmp -no-ansi-alias -restrict -xCORE-AVX2 -O3 -no-prec-div -ipo -inline-calloc -qopt-calloc
CFLAGS := -O3
CFLAGS := $(CFLAGS) -xCORE-AVX2
CFLAGS := $(CFLAGS) -no-prec-div
CFLAGS := $(CFLAGS) -qopenmp
CFLAGS := $(CFLAGS) -restrict
CFLAGS := $(CFLAGS) -no-ansi-alias
#CFLAGS := $(CFLAGS) -ansi-alias
CFLAGS := $(CFLAGS) -ipo
CFLAGS := $(CFLAGS) -inline-calloc
CFLAGS := $(CFLAGS) -qopt-calloc
CFLAGS := $(CFLAGS) -Wall
#CFLAGS := $(CFLAGS) -qopenmp-link=static
#CFLAGS := $(CFLAGS) -static-intel -wd10237
#CFLAGS := $(CFLAGS) -prof-gen
#CFLAGS := $(CFLAGS) -prof-use
#CFLAGS := $(CFLAGS) -no-vec -no-simd
#CFLAGS := $(CFLAGS) -qopt-report=4 -qopt-report-phase=vec
#CFLAGS := $(CFLAGS) -w3
#CFLAGS := $(CFLAGS) -debug
CFLAGS := $(CFLAGS) -wd3948   # suppress warning about "pragma simd" deprecation

all: mbir_ct Gen_SysMatrix clean

clean:
	rm *.o 

OBJ = initialize.o recon3d.o heap.o icd3d.o A_comp.o allocate.o MBIRModularUtils.o

mbir_ct: mbir_ct.o $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@ -lm
	mv $@ ../bin

Gen_SysMatrix: Gen_SysMatrix.o $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@ -lm
	mv $@ ../bin

#%.o: %.c 
#	$(CC) -c $(CFLAGS) $< -o $@
