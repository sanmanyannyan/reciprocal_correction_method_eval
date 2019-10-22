CC=gcc
CFLAGS=-O3 -march=native -mtune=cortex-a72
LDFLAGS=-lm

run:div
	./div > out.txt
div:div.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)
div.s:div.c
	$(CC) -S $(CFLAGS) $< -o $@ $(LDFLAGS)

.PHONY:clean run
clean:
	rm -f div.s div *.txt
