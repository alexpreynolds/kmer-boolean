SHELL := /bin/bash
PWD := $(shell pwd)
CC = gcc
CXX = g++
UNAME := $(shell uname -s)
BLDFLAGS = -Wall -Wextra -std=c++14
BLDDFLAGS = -Wall -Wextra -std=c++14 -pedantic
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
	
debug:
	$(CXX) -g $(BLDDFLAGS) $(CXXDFLAGS) -c kmer-boolean.cpp -o kmer-boolean.o
	$(CXX) -g $(BLDDFLAGS) $(CXXDFLAGS) -I$(INCLUDES) kmer-boolean.o -o kmer-boolean
		
test: cpp_test

cpp_test: cpp_test_A cpp_test_B cpp_test_C

cpp_test_A: cpp_test_A_sequential

cpp_test_A_sequential:
	${PWD}/kmer-boolean --all --k=2 --read-in-all-at-once < test/testA.fa > test/testA.obs
	@diff -s test/testA.obs test/testA.exp
	@rm test/testA.obs
	
cpp_test_B: cpp_test_B_sequential

cpp_test_B_sequential:
	${PWD}/kmer-boolean --all --k=2 --read-in-all-at-once < test/testB.fa > test/testB.obs
	@diff -s test/testB.obs test/testB.exp
	@rm test/testB.obs
	
cpp_test_C: cpp_test_C_stream

cpp_test_C_stream:
	${PWD}/kmer-boolean --present --k=2 < test/testC.fa > test/testC.present.obs
	${PWD}/kmer-boolean --absent --k=2 < test/testC.fa > test/testC.absent.obs
	${PWD}/kmer-boolean --all --k=2 < test/testC.fa > test/testC.all.obs
	sort <(cat test/testC.present.obs test/testC.absent.obs) > test/testC.obs	
	@diff -s test/testC.obs test/testC.exp
	@diff -s <(sort test/testC.all.obs) test/testC.exp
	@rm test/testC.obs test/testC.present.obs test/testC.absent.obs test/testC.all.obs
	
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