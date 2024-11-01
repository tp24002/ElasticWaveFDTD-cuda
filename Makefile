test.out:	test.o init.o update.o parameter.o memory.o
		nvcc -rdc=true -o test.out ./object/test.o ./object/init.o ./object/update.o ./object/memory.o ./object/parameter.o -Xcompiler -fopenmp -O3
# sample.out:	main.o init.o update.o parameter.o
# 		nvcc -rdc=true -o sample.out ./object/main.o ./object/init.o ./object/update.o ./object/parameter.o -Xcompiler -fopenmp -O3
# main.o:	./main.cu
# 		nvcc -rdc=true -o ./object/main.o -c ./main.cu
test.o:	./test1.cu
		nvcc -rdc=true -o ./object/test.o -c ./test1.cu
memory.o:	./src/memory.cu
		nvcc -o ./object/memory.o -c ./src/memory.cu
init.o: ./src/init.cu
		nvcc -rdc=true -o ./object/init.o -c ./src/init.cu
update.o: ./src/update.cu
		nvcc -rdc=true -o ./object/update.o -c ./src/update.cu -Xcompiler -fopenmp -O3
parameter.o: ./src/parameter.cu
		nvcc -rdc=true -o ./object/parameter.o -c ./src/parameter.cu