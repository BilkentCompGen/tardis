CC=gcc
CFLAGS = -O3 -g -I ../htslib -mcmodel=medium
LDFLAGS = ../htslib/libhts.a  -lz -lm -lpthread
SOURCES = vh_conflict.c vh_buffer.c vh_setcover.c vh_main.c vh_logger.c vh_divethandler.c vh_hash.c vh_createMaxClusterDeletion.c vh_heap.c vh_intervalhandler.c vh_maximalCluster.c vh_createMaxClusterInversion.c vh_createMaxClusterTDup.c vh_createMaxClusterMEI.c vh_createMaxClusterInsertion.c
	
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = libvh.a
INSTALLPATH = /usr/local/bin/

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS) 
	ar -rc libvh.a  $(OBJECTS) 

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean: 
	rm -f $(EXECUTABLE) *.o *~ 
