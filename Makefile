TARDIS_VERSION := "0.2 - The Universe Reboot Edition"
TARDIS_UPDATE := "June 26, 2017"
TARDIS_DEBUG := 1
BUILD_DATE := "$(shell date)"
CC=gcc
CFLAGS =  -O3 -g -I htslib -I vh -I sonic -DTARDIS_VERSION=\"$(TARDIS_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DTARDIS_UPDATE=\"$(TARDIS_UPDATE)\" -DTARDIS_DEBUG=$(TARDIS_DEBUG)
LDFLAGS = htslib/libhts.a vh/libvh.a sonic/libsonic.a -lz -lm -lpthread
SOURCES = tardis.c cmdline.c common.c processbam.c config.c processfq.c external.c vh.c variants.c splitread.c refgenome.c bamonly.c 
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
	rm -f $(EXECUTABLE) *.o *~

libs:
	make -C htslib
	make -C vh
	make -C sonic

install:
	cp tardis $(INSTALLPATH)
