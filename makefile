GCC=g++
OPTS=-O0

default:
	$(GCC) tests.cpp -o test.out $(OPTS)
	./test.out
	$(GCC) main.cpp -o main.out $(OPTS)

verbose:
	$(GCC) tests.cpp -o test.out -DVERBOSE $(OPTS)
	./test.out
	$(GCC) main.cpp -o main.out -DVERBOSE $(OPTS)
notest:
	$(GCC) main.cpp -o main.out $(OPTS)

debug:
	$(GCC) tests.cpp -g -o test.out $(OPTS) -DDEBUG
	./test.out
	$(GCC) main.cpp -g -o main.out $(OPTS) -DDEBUG

notestdebug:
	$(GCC) main.cpp -g -o main.out $(OPTS) -DDEBUG
	
