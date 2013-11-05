GCC=gcc
OPTS=-lm -O0
SOURCEF=main.c
OUT=-o main.out

default:
	$(GCC) $(SOURCEF) $(OUT) $(OPTS)


debug:
	$(GCC) $(SOURCEF) $(OUT) -g $(OPTS) -DDEBUG

snow:
	$(GCC) $(SOURCEF) $(OUT) $(OPTS) -DWSNOW

snowdebug:
	$(GCC) $(SOURCEF) $(OUT) $(OPTS) -DWSNOW -g -DDEBUG
	
fast:
	$(GCC) $(SOURCEF) -o fast.out -ffast-math -O3
	
