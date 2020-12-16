main: main.o model.o
	g++ -O3 -fopenmp -o meiosis main.o model.o -g

main.o: main.cpp model.h
	g++ -O3 -fopenmp -o main.o -c main.cpp -g

model.o: model.cpp model.h
	g++ -O3 -fopenmp -o model.o -c model.cpp -g

clean:
	rm main main.o model.o
