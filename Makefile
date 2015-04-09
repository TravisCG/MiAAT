all:
	gcc miaat.c -Wall -W -pedantic -lm -o miaat
clean:
	rm miaat *~
