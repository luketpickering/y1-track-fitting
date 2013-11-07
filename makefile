GCC=gcc
OPTS=-std=gnu99 -lc -lm
SOURCEF=main.c
EXE=main.out
OUT=-o $(EXE)
RTARGSF=args.txt
RUNOUTF=out.dat
EXEREDIR=>$(RUNOUTF)

default: defargs custargs standard
	
custargs: min_py standard
	./$(EXE) `cat $(RTARGSF)` $(EXEREDIR)
	echo 'IAMSILLY=True' > settings.py
	./hist.py $(RUNOUTF)

defargs: $(RTFARGSF)
	echo '-f manytracks.raw -m 5 -s 0.05' > $(RTARGSF)

min_py: __init__.py settings.py
	touch __init__.py
	touch settings.py
	
expl:
	pdflatex ./expl-src/expl.tex

standard:
	$(GCC) $(SOURCEF) $(OUT) $(OPTS) -O2 -DWSNOW

nosnow:
	$(GCC) $(SOURCEF) $(OUT) $(OPTS) -O2
fast:
	$(GCC) $(SOURCEF) -o fast.out -ffast-math -O3

debug:
	$(GCC) $(SOURCEF) $(OUT) -g $(OPTS) -DDEBUG

notsilly: min_py
	$(GCC) $(SOURCEF) $(OUT) $(OPTS) -O2
	echo 'IAMSILLY=False' > settings.py
	

	

	
