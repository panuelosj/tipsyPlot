CC = gcc
CFLAGS = -c -Wall
LDFLAGS =
LIBS = -lm -lX11 -lplplotd -lconfig

# files
SRD_SRC = $(wildcard shared/*.c)
SRD_OBJ = $(SRD_SRC:.c=.o)
EXEC = tipsyPlot
EXEC2 = plotSlice
EXEC3 = plotTempSlice
EXEC4 = tipsy2gadget
EXEC5 = gadgetPlot
EXEC6 = gadgetSlice
EXEC7 = gadgetTempSlice
EXEC8 = plotRhoSlice
EXEC9 = plotDumpSlice

# main
$(EXEC): $(EXEC).o $(SRD_OBJ) $(EXEC).h
	$(CC) $(LDFLAGS) $? $(LIBS) -o $@

$(EXEC2): $(EXEC2).o $(SRD_OBJ) $(EXEC).h
	$(CC) $(LDFLAGS) $? $(LIBS) -o $@

$(EXEC3): $(EXEC3).o $(SRD_OBJ) $(EXEC).h
	$(CC) $(LDFLAGS) $? $(LIBS) -o $@

$(EXEC4): $(EXEC4).o $(SRD_OBJ) $(EXEC).h
	$(CC) $(LDFLAGS) $? $(LIBS) -o $@

$(EXEC5): $(EXEC5).o $(SRD_OBJ) $(EXEC).h
	$(CC) $(LDFLAGS) $? $(LIBS) -o $@

$(EXEC6): $(EXEC6).o $(SRD_OBJ) $(EXEC).h
	$(CC) $(LDFLAGS) $? $(LIBS) -o $@

$(EXEC7): $(EXEC7).o $(SRD_OBJ) $(EXEC).h
	$(CC) $(LDFLAGS) $? $(LIBS) -o $@

$(EXEC8): $(EXEC8).o $(SRD_OBJ) $(EXEC).h
	$(CC) $(LDFLAGS) $? $(LIBS) -o $@

$(EXEC9): $(EXEC9).o $(SRD_OBJ) $(EXEC).h
	$(CC) $(LDFLAGS) $? $(LIBS) -o $@

%.o: %.c $(EXEC).h
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -v *.o
	rm -v shared/*.o
