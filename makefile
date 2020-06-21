CC=c++
STD=--std=c++17

make:main.cpp
	$(CC) -I -w $(STD) -o main.o main.cpp

clean:
	rm main
