SHELL=/bin/sh
OS = MACOSX


COMPILER=g++
OPT=-Wall -O2 -D$(OS)

HEADER= locus.h sample.h population.h random.h coalescent_tree.h coalescent_tree_diversity.h distribution.h readfile.h matrix.h


#infoseq.h infoseq_diversity.h 

SRC= locus.cpp\
     sample.cpp\
     coalescent_tree.cpp\
     coalescent_tree_diversity.cpp\
     population.cpp\
     random.cpp\
     distribution.cpp\
     distributions.cpp\
     readfile.cpp\
     locus_TE.cpp

SRC_C= matrix.c
 

OBJ= $(SRC:.cpp=.o) $(SRC_C:.c=.o)

BIN= Tests_allcoal\
     SiteFrequencySpectrum\
     Tests_variance\
     SpecCoalABGD\
     CoalABGD\
     CI_Theta_estimators\
     CoalTE\
     generate_CI_coal\
     GenerCoalSeq\
     SkylinePlot\
     SerialCoalSeq


all: $(BIN)

test: main_test.cpp $(OBJ) $(HEADER)
	$(COMPILER) $(OPT) -o $@ $(OBJ) main_test.cpp -lm

tst2: main_tst2.cpp $(OBJ) $(HEADER)
	$(COMPILER) $(OPT) -o $@ $(OBJ) main_tst2.cpp -lm

SkewedTree: main_SkewedTree.cpp $(OBJ) $(HEADER)
	$(COMPILER) $(OPT) -o $@ $(OBJ) main_SkewedTree.cpp -lm

CoalTst: main_tst.cpp $(OBJ) $(HEADER)
	$(COMPILER) $(OPT) -o $@ $(OBJ) main_tst.cpp -lm

CoalTE: main_CoalTE.cpp $(OBJ) $(HEADER)
	$(COMPILER) $(OPT) -o $@ $(OBJ) main_CoalTE.cpp -lm

Tests_allcoal: main_Tests_All_Coal.cpp $(OBJ) $(HEADER)
	$(COMPILER) $(OPT) -o $@ $(OBJ) main_Tests_All_Coal.cpp -lm

SiteFrequencySpectrum: main_SiteFrequencySpectrum.cpp $(OBJ) $(HEADER)
	$(COMPILER) $(OPT) -o $@ $(OBJ) main_SiteFrequencySpectrum.cpp -lm

Tests_variance: main_Tests_Var.cpp $(OBJ) $(HEADER)
	$(COMPILER) $(OPT) -o $@ $(OBJ) main_Tests_Var.cpp -lm
	
SpecCoal: main_SpecCoal.cpp  $(OBJ) $(HEADER)
	$(COMPILER) $(OPT) -o $@ $(OBJ) main_SpecCoal.cpp

SpecCoalABGD: main_SpecCoalABGD.cpp abgd.c $(OBJ) $(HEADER)
	$(COMPILER) $(OPT) -o $@ $(OBJ) main_SpecCoalABGD.cpp

CoalABGD: main_CoalABGD.cpp $(OBJ) $(HEADER)
	$(COMPILER) $(OPT) -o $@ $(OBJ) main_CoalABGD.cpp

CI_Theta_estimators: main_CI_Theta_estimators.cpp $(OBJ) $(HEADER)
	$(COMPILER) $(OPT) -o $@ $(OBJ) main_CI_Theta_estimators.cpp

generate_CI_coal: main_generate_CI_coal.cpp $(OBJ) $(HEADER)
	$(COMPILER) $(OPT) -o $@ $(OBJ) main_generate_CI_coal.cpp

GenerCoalSeq: main_GenerCoalSeq.cpp $(OBJ) $(HEADER)
	$(COMPILER) $(OPT) -o $@ $(OBJ) main_GenerCoalSeq.cpp

CoalDiv: main_CoalDiv.cpp $(OBJ) $(HEADER)
	$(COMPILER) $(OPT) -o $@ $(OBJ) main_CoalDiv.cpp

CoalSplit: main_CoalSplit.cpp $(OBJ) $(HEADER)
	$(COMPILER) $(OPT) -o $@ $(OBJ) main_CoalSplit.cpp

SkylinePlot: main_SkylinePlot.cpp $(OBJ) $(HEADER)
	$(COMPILER) $(OPT) -o $@ $(OBJ) main_SkylinePlot.cpp

SerialCoalSeq: main_SerialCoalSeq.cpp $(OBJ) $(HEADER)
	$(COMPILER) $(OPT) -o $@ $(OBJ) main_SerialCoalSeq.cpp


%.o: %.cpp
	$(COMPILER) $(OPT) -c $<

%.o: %.cpp %.h
	$(COMPILER) $(OPT) -c $<

%.o: %.c %.h
	$(COMPILER) $(OPT) -c $<


clean:
	\rm -f *.o $(BIN)


binclean:
	cd ~/bin
	\rm -f $(BIN)


install: all
	cp $(BIN) ~/bin


archive: clean
	\rm -rf SimulTrees; mkdir SimulTrees;
	cp -r *.cpp *.h *.c Makefile SimulTrees ;
	tmp=`date | awk '{print $$3$$2$$NF}'`; tar czvf ../SimulTrees.$$tmp.tgz SimulTrees;
	\rm -rf SimulTrees

