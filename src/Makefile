
.PHONY:all clean install 

CC=gcc
CFLAGS=-O3 
CFLAGS+=-Wall -Wno-unused-result -Wno-unused-function -Werror -Werror=unused-but-set-variable
LDLIBS=-lz
ALL=shmr_mkseqdb shmr_index shmr_gather_mc shmr_overlap shmr_dedup shmr_map shmr_overlap_annotate

all: $(ALL)

kalloc.o: kalloc.c

shmr_reduce.o: shmr_reduce.c

mm_sketch.o: mm_sketch.c

shmr_utils.o: shmr_utils.c

shmr_utils.so: shmr_utils.c kalloc.c
	gcc -O3 -shared -fPIC -Wall kalloc.c shmr_utils.c -o shmr_utils.so

shimmer4py.so: shimmer4py.c shmr_utils.c kalloc.c
	gcc -O3 -shared -fPIC -Wall kalloc.c shmr_utils.c shimmer4py.c -o shimmer4py.so

shmr_end_filter.o: shmr_end_filter.c

shmr_gather_mc.o: shmr_gather_mc.c

DWmatch.o: DWmatch.c

shmr_mkseqdb: shmr_mkseqdb.c kalloc.o shmr_utils.o

shmr_index: shmr_index.c kalloc.o shmr_reduce.o mm_sketch.o shmr_utils.o shmr_end_filter.o

shmr_gather_mc: shmr_gather_mc.o kalloc.o shmr_utils.o

shmr_overlap: shmr_overlap.c shmr_utils.o kalloc.o  DWmatch.o

shmr_map: shmr_map.c shmr_utils.o kalloc.o  DWmatch.o

shmr_dedup: shmr_dedup.c kalloc.o

shmr_overlap_annotate: shmr_overlap_annotate.c kalloc.o shmr_utils.o

shmr_dedup2: shmr_dedup2.c kalloc.o

BINDIR=$(shell dirname $(shell which python))
install:
	cp $(ALL) $(BINDIR) 

clean:
	rm -f shmr_dedup shmr_gather_mc shmr_overlap \
		shmr_map shmr_mkseqdb shmr_index shmr_overlap_annotate \
		shmr_gather_mc *.o *.so; rm -rf ./bin/
