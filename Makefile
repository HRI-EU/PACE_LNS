CC=g++
CFLAGS = -O3 -std=c++17 -Wno-write-strings -Wno-unused-result
LDFLAGS =
INCLUDES =
TARGET = main

main: main.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCLUDES) main.cpp -o main


clean:
	rm -f main
	rm -rf *.o
