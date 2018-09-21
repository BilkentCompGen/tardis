TARDIS_VERSION := "0.5.1"
TARDIS_UPDATE := "August 29, 2018"
TARDIS_DEBUG := 0
BUILD_DATE := "$(shell date)"
CC=gcc
CFLAGS =  -O3 -g -I htslib -I vh -I sonic -DTARDIS_VERSION=\"$(TARDIS_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DTARDIS_UPDATE=\"$(TARDIS_UPDATE)\" -DTARDIS_DEBUG=$(TARDIS_DEBUG) 
LDFLAGS = htslib/libhts.a vh/libvh.a sonic/libsonic.a -lz -lm -lpthread -llzma -lbz2 
SOURCES = tardis.c cmdline.c common.c processbam.c config.c processfq.c external.c variants.c splitread.c processrefgen.c bamonly.c free.c mappings.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = tardis
INSTALLPATH = /usr/local/bin/

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	make clean -C htslib
	make clean -C sonic
	make clean -C vh
	rm -f $(EXECUTABLE) *.o *~

libs:
	make -C htslib
	make -C vh
	make -C sonic

install:
	cp tardis $(INSTALLPATH)
