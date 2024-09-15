sample.out:	main.o init.o insert.o update.o para_in.o extradition.o
		nvcc -o sample.out ./object/main.o ./object/init.o ./object/insert.o ./object/update.o ./object/extradition.o ./object/para_in.o -Xcompiler -fopenmp -O3
# sample.out:	main2.o init.o insert.o update.o print.o
# 			gcc -o sample.out ./object/main2.o ./object/init.o ./object/insert.o ./object/update.o ./object/print.o -lm -fopenmp -O3
main.o:	./main.cu
		nvcc -o ./object/main.o -c ./main.cu
extradition.o:	./src/extradition.cu
		nvcc -o ./object/extradition.o -c ./src/extradition.cu
init.o: ./src/init.cu
		nvcc -o ./object/init.o -c ./src/init.cu
insert.o: ./src/insert.cu
		nvcc -o ./object/insert.o -c ./src/insert.cu -lm
update.o: ./src/update.cu
		nvcc -o ./object/update.o -c ./src/update.cu -Xcompiler -fopenmp -O3
# print.o: ./src/print.cu
# 		nvcc -o ./object/print.o -c ./src/print.cu
para_in.o: ./src/para_in.cu
		nvcc -o ./object/para_in.o -c ./src/para_in.cu