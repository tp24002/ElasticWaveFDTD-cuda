sample.out:	main.o init.o insert.o update.o para_in.o print.o 
		gcc -o sample.out ./object/main.o ./object/init.o ./object/insert.o ./object/update.o ./object/print.o ./object/para_in.o -lm -fopenmp -O3
# sample.out:	main2.o init.o insert.o update.o print.o
# 			gcc -o sample.out ./object/main2.o ./object/init.o ./object/insert.o ./object/update.o ./object/print.o -lm -fopenmp -O3
main.o:	./main.c
		gcc -o ./object/main.o -c ./main.c
# main.o:	./main2.c
# 		gcc -o ./object/main.o -c ./main2.c
# main.o: ./main3.c
# 		gcc -o ./object/main.o -c ./main3.c
init.o: ./src/init.c
		gcc -o ./object/init.o -c ./src/init.c
insert.o: ./src/insert.c
		gcc -o ./object/insert.o -c ./src/insert.c -lm
update.o: ./src/update.c
		gcc -o ./object/update.o -c ./src/update.c -lm -fopenmp -O3
print.o: ./src/print.c
		gcc -o ./object/print.o -c ./src/print.c
para_in.o: ./src/para_in.c
		gcc -o ./object/para_in.o -c ./src/para_in.c