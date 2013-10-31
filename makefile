GCC=g++
OPTS=-O0 -DWSNOW

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

notestfast:
	$(GCC) main.cpp -o fast.out -O1
debug:
	$(GCC) tests.cpp -g -o test.out $(OPTS) -DDEBUG
	./test.out
	$(GCC) main.cpp -g -o main.out $(OPTS) -DDEBUG

notestdebug:
	$(GCC) main.cpp -g -o main.out $(OPTS) -DDEBUG
	
