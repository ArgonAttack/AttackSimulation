CC = gcc
CFLAGS = -Wall -Wextra -Wno-unused-variable -Wno-unused-parameter
CFLAGSRUN = -O2
CFLAGSDEBUG = -O0 -g
LDFLAGS = -lm
SOURCES = Attack.c BuildGraph.c statistics.c main.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE = Attack


all: $(EXECUTABLE)
    
Attack: $(SOURCES)
	$(CC) $(CFLAGS) $(CFLAGSRUN) -o $(EXECUTABLE) $(SOURCES) $(LDFLAGS)

debug: $(SOURCES)
	$(CC) $(CFLAGS) $(CFLAGSDEBUG) -o $(EXECUTABLE) $(SOURCES) $(LDFLAGS)

clean:
	rm Attack