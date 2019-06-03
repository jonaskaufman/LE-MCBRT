CC = nvcc

# load cuda in the shell prompt
#	module load cuda
	
all: ray.o simulation.cu
	$(CC) -std=c++11 ray.o -o simulation.x simulation.cu

ray.o: ray.cpp
	$(CC) -std=c++11 -c ray.cpp

clean:
	rm -f *.o *.x
