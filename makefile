SHELL := /bin/bash
PWD := $(shell pwd)
CC = gcc
CXX = g++
UNAME := $(shell uname -s)
BLDFLAGS = -Wall -Wextra -std=c++14
BLDDFLAGS = -Wall -Wextra -std=c++14 -pendantic
CXXFLAGS = -D__STDC_CONSTANT_MACROS -D__STDINT_MACROS -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE=1 -O3
CXXDFLAGS = -D__STDC_CONSTANT_MACROS -D__STDINT_MACROS -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE=1 -O3 -DDEBUG_FLAG=1
INCLUDES = /usr/include

ifeq ($(UNAME),Darwin)
	CC = clang
	CXX = clang++
	FLAGS += -Weverything
else
	CXXFLAGS += -static
endif

all: kmer-boolean

kmer-boolean:
	$(CXX) -g $(BLDFLAGS) $(CXXFLAGS) -c kmer-boolean.cpp -o kmer-boolean.o
	$(CXX) -g $(BLDFLAGS) $(CXXFLAGS) -I$(INCLUDES) kmer-boolean.o -o kmer-boolean

clean:
	rm -rf *~
	rm -rf kmer-boolean
	rm -rf kmer-boolean.o

python: python_test

python_test: python_true python_false

python_true:
	echo -e ">foo\nCATTCTC\nGGGAC\n>bar\nTTATAT\n>baz\nTTTATTAG\nACCTCT" | ./kmer-boolean.py 2 CG

python_false:
	-echo -e ">foo\nCATTCTC\nGGGAC\n>bar\nTTATAT\n>baz\nTTTATTAG\nACCTCT" | ./kmer-boolean.py 2 AA