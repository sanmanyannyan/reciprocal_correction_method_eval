CC=gcc
CFLAGS=-O3 -mfpu=neon -march=native
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
