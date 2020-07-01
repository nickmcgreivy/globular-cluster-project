CC=g++

make:main.cpp
	@$(CC) -w main.cpp
	@./a.out
clean:
	@rm a.out
