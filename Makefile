CC = gcc -Wall
CFLAGS = -Wall -O3 -std=c99 -lm
# Add -g for debugging, remove -O3 for faster compilation during development
# CFLAGS = -Wall -g -std=c99 -lm

SRCDIR = src
INCLUDEDIR = include
BUILDDIR = build
BINDIR = $(BUILDDIR)/bin

SOURCES = $(SRCDIR)/main.c \
          $(SRCDIR)/calc.c \
          $(SRCDIR)/io.c \
          $(SRCDIR)/boundaries.c \
          $(SRCDIR)/memory.c \
          $(SRCDIR)/velocity.c

HEADERS = $(INCLUDEDIR)/constants.h \
          $(INCLUDEDIR)/variables.h \
          $(INCLUDEDIR)/initialise.h \
          $(INCLUDEDIR)/functions.h \
          $(INCLUDEDIR)/calc.h \
          $(INCLUDEDIR)/velocity.h \
          $(INCLUDEDIR)/memory.h

OBJECTS = $(patsubst $(SRCDIR)/%.c,$(BUILDDIR)/%.o,$(SOURCES))

TARGET = $(BINDIR)/srd_simulation

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $^ -o $@

$(BUILDDIR)/%.o: $(SRCDIR)/%.c $(HEADERS)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -I$(INCLUDEDIR) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	@rm -rf $(BUILDDIR)
	@rm -f data/*.dat data/*.pdb
