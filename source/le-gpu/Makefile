# load cuda in the shell prompt
#	module load cuda 

objects = simulation.o ray.o

all: $(objects)
	nvcc -arch=sm_30 -Xptxas -O1 $(objects) -g -o simulation.x

%.o : %.cu
	nvcc -std=c++11 -Xptxas -O1 -lcurand -x cu -arch=sm_30 -I. -dc $< -o $@	

clean:
	rm -f *.o *.x *.csv
