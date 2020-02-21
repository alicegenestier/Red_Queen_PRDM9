main: main.o model.o
	g++ -o meiosis main.o model.o -g

main.o: main.cpp model.h
	g++ -o main.o -c main.cpp -g

model.o: model.cpp model.h
	g++ -o model.o -c model.cpp -g

clean:
	rm main main.o model.o
