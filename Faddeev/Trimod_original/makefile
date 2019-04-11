FC = gfortran

OPT = -O3 -llapack -lblas

RUNO    = gauleg.o gauss.o spline.o pot.o tmat.o search.o main.o \

all: run

run: $(RUNO)
	${FC} ${OPT} -o run $(RUNO)

gauleg.o            : gauleg.f
					${FC} -c gauleg.f

gauss.o             : gauss.f
					${FC} -c gauss.f

spline.o            : spline.f
					${FC} -c spline.f

pot.o               : pot.f
					${FC} -c pot.f

tmat.o              : tmat.f
					${FC} -c tmat.f

search.o            : search.f
					${FC} -c search.f

main.o              : main.f
					${FC} -c main.f

# For a fresh build 
clean:
	rm -f *.o *.mod *.MOD
