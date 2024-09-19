sample.out:	main.o init.o insert.o update.o parameter.o memory.o
		nvcc -o sample.out ./object/main.o ./object/init.o ./object/insert.o ./object/update.o ./object/memory.o ./object/parameter.o -Xcompiler -fopenmp -O3
# sample.out:	main2.o init.o insert.o update.o print.o
# 			gcc -o sample.out ./object/main2.o ./object/init.o ./object/insert.o ./object/update.o ./object/print.o -lm -fopenmp -O3
main.o:	./main.cu
		nvcc -o ./object/main.o -c ./main.cu
memory.o:	./src/memory.cu
		nvcc -o ./object/memory.o -c ./src/memory.cu
init.o: ./src/init.cu
		nvcc -o ./object/init.o -c ./src/init.cu
insert.o: ./src/insert.cu
		nvcc -o ./object/insert.o -c ./src/insert.cu -lm
update.o: ./src/update.cu
		nvcc -o ./object/update.o -c ./src/update.cu -Xcompiler -fopenmp -O3
# print.o: ./src/print.cu
# 		nvcc -o ./object/print.o -c ./src/print.cu
parameter.o: ./src/parameter.cu
		nvcc -o ./object/parameter.o -c ./src/parameter.cu