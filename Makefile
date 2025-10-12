# компиляция под АРМ
# export CROSS_COMPILE=arm-linux-gnueabihf-
# PKG_CONFIG_PATH=C:\GTK\lib/pkgconfig
# ARCH ?= arm
# CROSS_COMPILE ?= arm-none-eabi-


ifeq ($(ARCH),arm)
OPTIMIZATION = -mtune=cortex-a8 -mfpu=neon -mfloat-abi=hard
#OPTIMIZATION = -march=armv7-a -mtune=cortex-a8 -mfpu=vfpv3 -mfloat-abi=hard
else
OPTIMIZATION += -march=native

ARCH ?= $(shell uname -m)
ifeq ($(ARCH),aarch64)
else ifeq ($(ARCH),x86_64)
else
endif
endif
OPTIMIZATION += -ffast-math -ftree-vectorize

CXX= $(CROSS_COMPILE)g++ $(OPTIMIZATION)
CC = $(CROSS_COMPILE)gcc $(OPTIMIZATION)
AS = $(CROSS_COMPILE)gcc $(OPTIMIZATION)
SIZE = $(CROSS_COMPILE)size
OBJCOPY = $(CROSS_COMPILE)objcopy
CFLAGS = -Werror -std=gnu11 -D_GNU_SOURCE
CXXFLAGS = -Werror -std=c++17 -D_GNU_SOURCE
CFLAGS += -O3 
# $(shell pkg-config --cflags )
CXXFLAGS += -O3
LDFLAGS = -g
# -g $(shell pkg-config --libs   ) -lOpenCL
# Идентификация операционной системы
UNAME=$(shell uname -o)

# в зависимости от операционной системы
ifeq ($(UNAME), Msys)
else ifeq ($(UNAME), GNU/Linux)
endif
# SRC_XX = kernel.cpp main.cpp 

SRC = lzjb.c lzjb_original.c huffman.c deflate.c gzip.c xxh64.c crc.c

OUTPUT=lzjb2

ASM_OBJECTS := $(ASMSRC:.s=.o)


all: $(OUTPUT)

#	wc -l $(SRC)
$(OUTPUT): $(ASM_OBJECTS) $(SRC:.c=.o) $(SRC_XX:.cpp=.o)
	$(CXX) -o $(OUTPUT).exe $^ $(LDFLAGS)
	wc -l $(SRC) $(SRC_XX)
	$(SIZE) $^ $(OUTPUT).exe

$(ASM_OBJECTS): %.o : %.s
	$(CC) $(ASFLAGS) -c -o $@ $<

$(SRC_XX:.cpp=.o): %.o : %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OUTPUT) $(SRC:.c=.o) $(SRC_XX:.cpp=.o) $(ASM_OBJECTS)

