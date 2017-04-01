#
# Makefile for smooth.
#
CFLAGS	=   -g
LIBS	=   -lm

default:	spherematchDJE 

clean:
	rm -f *.o

spherematchDJE: main.o kd.o smooth.o
	$(CC) $(CFLAGS) -o spherematchDJE main.o kd.o smooth.o $(LIBS)

main.o: main.c kd.h smooth.h 
	$(CC) $(CFLAGS) -c main.c

kd.o: kd.c kd.h 
	$(CC) $(CFLAGS) -c kd.c

smooth.o: smooth.c kd.h smooth.h
	$(CC) $(CFLAGS) -c smooth.c

tar:
	tar cvf spherematchDJE.tar Makefile kd.c kd.h main.c smooth.c smooth.h 
