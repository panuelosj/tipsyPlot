CC = gcc
CFLAGS = -c -Wall
LDFLAGS =
LIBS = -lm -lX11 -lplplotd

# files
SRD_SRC = $(wildcard shared/*.c)
SRD_OBJ = $(SRD_SRC:.c=.o)
EXEC = tipsyPlot

# main
$(EXEC): $(EXEC).o $(SRD_OBJ) $(EXEC).h
	$(CC) $(LDFLAGS) $? $(LIBS) -o $@

%.o: %.c $(EXEC).h
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -v *.o
	rm -v shared/*.o
