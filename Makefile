main: main.o model.o
	g++ -o meiosis main.o model.o

main.o: main.cpp model.h
	g++ -o main.o -c main.cpp

model.o: model.cpp model.h
	g++ -o model.o -c model.cpp

clean:
	rm main main.o model.o
