# load cuda in the shell prompt
# module load cuda 

EXEC = base-run.x le-run.x

all: $(EXEC)

%-run.x: %-simulation.o ray.o
	nvcc -arch=sm_30 $^ -g -o $@

%.o : %.cu
	nvcc -arch=sm_30 -std=c++11 -lcurand -x cu -I. -dc $< -o $@	

clean:
	rm -f *.o *.x *.csv
