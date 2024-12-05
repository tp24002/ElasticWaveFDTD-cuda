# test.out:	test.o init.o update.o parameter.o memory.o
# 		nvcc -rdc=true -o test.out ./object/test.o ./object/init.o ./object/update.o ./object/memory.o ./object/parameter.o -Xcompiler -fopenmp -O3

sample.out:	main.o init.o update.o parameter.o memory.o visualization.o
		nvcc -rdc=true -o sample.out ./object/main.o ./object/init.o ./object/update.o ./object/memory.o ./object/parameter.o ./object/visualization.o -Xcompiler -fopenmp -O3 -lpng
main.o:	./main.cu
		nvcc -rdc=true -o ./object/main.o -c ./main.cu
init.o: ./src/init.cu
		nvcc -rdc=true -o ./object/init.o -c ./src/init.cu
update.o: ./src/update.cu
		nvcc -rdc=true -o ./object/update.o -c ./src/update.cu -Xcompiler -fopenmp -O3
parameter.o: ./src/parameter.cu
		nvcc -rdc=true -o ./object/parameter.o -c ./src/parameter.cu
memory.o:	./src/memory.cu
		nvcc -rdc=true -o ./object/memory.o -c ./src/memory.cu
visualization.o:	./src/visualization.cu
		nvcc -x cu -o ./object/visualization.o -c ./src/visualization.cu -lpng
# test.o:	./test1.cu
# 		nvcc -rdc=true -o ./object/test.o -c ./test1.cu
