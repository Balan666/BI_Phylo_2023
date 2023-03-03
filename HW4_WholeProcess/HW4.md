## 1. Какие программы использовали для анализа? версии программ.

1. trimAl (http://trimal.cgenomics.org/downloads или https://anaconda.org/bioconda/trimal; если совсем не взлетает,
то https://ngphylogeny.fr/tools/tool/284/form).
2. RAxML-NG, https://github.com/amkozlov/raxml-ng/releases/tag/1.1.0, https://raxml-ng.vital-it.ch/#/
3. ModelTest-NG, https://github.com/ddarriba/modeltest/releases/tag/v0.1.7
RAxML
и ModelTest есть на https://www.phylo.org/index.php
4. IQ-Tree (http://www.iqtree.org/, если не взлетает, ищите внизу
странички web service или https://anaconda.org/bioconda/iqtree). Для практикума
№4 подойдёт любая версия, для практикума
№5 лучше iq-tree v2, поэтому имеет смысл ставить
самую новую.
5. (опция для выполнения бонусного
задания!) PartitionFinder: http://www.robertlanfear.com/partitionfinder/, https://github.com/brettc/partitionfinder/releases/tag/v2.1.1.
Очень подробная
инструкция по установке PartitionFinder на новые Ubuntu:
https://youtu.be/P77StHPIJ0k?t=14


## 2. Как вырезать плохие участки из выравнивания с помощью trimAl?
```
#trimal -in SUP35_aln_prank.best.fas -out SUP35_aln_prank.trim.fas
trimal -in SUP35_aln_prank.best.fas -out SUP35_aln_prank.trim.fas -automated1
```
## 3. Как подобрать модель эволюции в ModelTest (ModelTest-NG)?

```
modeltest-ng-static -i SUP35_aln_prank.trim.fas -o SUP35_aln_prank_trim_modeltest
```

## 4. Какая(ие) модель(и) эволюции была признана наиболее подходящей для нашего выравнивания?.
Давайте для единообразия ориентироваться на BIC.
```
--------------------------------------------------------------------------------

BIC       model              K            lnL          score          delta    weight
--------------------------------------------------------------------------------
       1  TIM3+G4            7     -8998.9374     18180.5947         0.0000    0.4140
       2  TrN+I+G4           7     -8999.7992     18182.3182         1.7235    0.1749
       3  TIM3+I+G4          8     -8996.0035     18182.3401         1.7454    0.1730
       4  TrN+G4             6     -9003.7976     18182.7017         2.1070    0.1444
       5  TIM2+I+G4          8     -8997.4609     18185.2549         4.6603    0.0403
       6  TIM2+G4            7     -9001.4442     18185.6083         5.0136    0.0338
       7  GTR+G4             9     -8995.4907     18188.9278         8.3331    0.0064
       8  TIM1+I+G4          8     -8999.5655     18189.4642         8.8695    0.0049
       9  GTR+I+G4          10     -8992.0690     18189.6979         9.1032    0.0044
      10  TIM1+G4            7     -9003.5780     18189.8758         9.2811    0.0040
--------------------------------------------------------------------------------
Best model according to BIC
---------------------------
Model:              TIM3+G4
lnL:                -8998.9374
Frequencies:        0.3401 0.1970 0.2260 0.2369
Subst. Rates:       1.7149 5.6649 1.0000 1.7149 14.5197 1.0000 
Inv. sites prop:    -
Gamma shape:        0.4034
Score:              18180.5947
Weight:             0.4140
---------------------------
```
## 5. Постройте ML-дерево в RAxML-NG, используя выбранную модель. (Приведите код или описание действий.)

```
raxml-ng --check --msa SUP35_aln_prank.trim.fas  --model TIM3+G4 --prefix SUP35_raxml_test
raxml-ng --msa SUP35_aln_prank.trim.fas --model TIM3+G4 --prefix SUP35_raxml --threads 2 --seed 222 
```
## 6. Отрисуйте полученное дерево (лучшее ML-дерево) и покажите рисунок.
+1 балл за скрипт на любом языке, который принимает на вход название файла и автоматически рисует качественное дерево (скрипт можно привести внутри отчёта, вряд ли он очень длинный).

```{r}
library(ggtree)
tr <- read.tree("SUP35_raxml.raxml.bestTree")
ggtree(tr) + geom_tiplab() + xlim(0,2)
```
```
#```{python}
#from Bio import AlignIO
#alignments = AlignIO.parse('SUP35_aln_prank.best.fas', "fasta")
#AlignIO.write(alignments, 'SUP35_aln.best.p.phy', 'phylip-relaxed')
#AlignIO.write(alignments, 'SUP35_aln.best.p.phy', 'phylip-sequential')
#```
```

## 7. Как выбрать модель в ModelFinder (можно через IQ-TREE)?

```
#iqtree -s SUP35_aln_prank.trim.fas -m MFP -pre SUP35_MF
iqtree2 -m MFP -s SUP35_aln_prank.trim.fas --prefix SUP35_MF2

iqtree2 -m TIM3+F+G4 -s SUP35_aln_prank.trim.fas --prefix SUP35_iqtree

iqtree2 -s SUP35_aln_prank.best.fas -pre SUP35_prank_unfilt
iqtree2 -s ../../Phylo-3-alignment/SUP35_10seqs_mafft.fa -pre SUP35_mafft
iqtree2 -s ../../Phylo-3-alignment/SUP35_10seqs_kalign.fa -pre SUP35_kalign
iqtree2 -s ../../Phylo-3-alignment/SUP35_10seqs_kalign.fa -pre SUP35_kalign -redo
iqtree2 -s ../../Phylo-3-alignment/SUP35_10seqs_clustalw.fa -pre SUP35_clustalw
```

## 8. Какая модель эволюции была признана наиболее подходящей для нашего выравнивания?
(Скопируйте нужную часть вывода программы).

```
Best-fit model according to BIC: TIM3+F+G4

List of models sorted by BIC scores: 

Model                  LogL         AIC      w-AIC        AICc     w-AICc         BIC      w-BIC
TIM3+F+G4         -8993.686   18035.372 -   0.0328   18035.972 -   0.0354   18170.092 +    0.732
TIM3+F+I+G4       -8991.335   18032.671 +    0.127   18033.321 +    0.133   18173.004 +    0.171
TN+F+G4           -9000.108   18046.216 - 0.000145   18046.768 -  0.00016   18175.323 +   0.0535
TN+F+I+G4         -8997.491   18042.982 -  0.00073   18043.582 - 0.000787   18177.702 -   0.0163
GTR+F+G4          -8990.601   18033.202 +    0.097   18033.905 +   0.0994   18179.149 -   0.0079
TIM2+F+G4         -8998.340   18044.680 - 0.000312   18045.280 - 0.000337   18179.399 -  0.00697
TIM3+F+I+I+R2     -8990.918   18033.836 +   0.0706   18034.539 +   0.0724   18179.783 -  0.00576
TIM2+F+I+G4       -8995.740   18041.480 -  0.00155   18042.131 -  0.00163   18181.814 -  0.00209
GTR+F+I+G4        -8988.267   18030.534 +    0.368   18031.291 +    0.367   18182.094 -  0.00181
TIM+F+G4          -8999.957   18047.915 - 6.19e-05   18048.515 - 6.68e-05   18182.635 -  0.00138
TN+F+I+I+R2       -8996.634   18043.269 - 0.000632   18043.919 - 0.000665   18183.602 - 0.000853
TIM+F+I+G4        -8997.340   18044.680 - 0.000312   18045.330 - 0.000329   18185.013 - 0.000421
TIM2+F+I+I+R2     -8994.830   18041.659 -  0.00141   18042.362 -  0.00145   18187.606 - 0.000115
```

Отличаются ли модели, выбранные ModelTest и ModelFinder, и насколько сильно?
Отличаются незначительно, в целом в выдаче есть много общего

## 9. Постройте ML-дерево в IQ-TREE, используя выбранную модель. (Приведите код или описание действий.)

в целом IQ-tree рисует его сразу в том же файле в ASCII:
```
+----------------------------------SUP35_Kla_AB039749
|
+-------------------------------------------SUP35_Agos_ATCC_10895_NM_211584
|
|                                                         +**SUP35_Scer_74-D694_GCA_001578265.1
|                                                      +**|
|                                                      |  +**SUP35_Scer_beer078_CM005938
|                                                +-----|
|                                                |     +**SUP35_Sbou_unique28_CM003560
|                                            +---|
|                                            |   +---SUP35_Spar_A12_Liti
|                                         +--|
|                                         |  +--------SUP35_Smik_IFO1815T_30
+-----------------------------------------|
                                          |     +---------SUP35_Sarb_H-6_chrXIII_CM001575
                                          |  +--|
                                          |  |  +------------SUP35_Seub_CBS12357_chr_II_IV_DF968535
                                          +--|
                                             +------------SUP35_Skud_IFO1802T_36
```

## 10. Отрисуйте полученное дерево (лучшее ML-дерево) и покажите рисунок.

```{r}
library(ggtree)
tr <- read.tree("/home/bananna/Desktop/BI_Phylo/HW4/SUP35_raxml.raxml.bestTree")
ggtree(tr) + geom_tiplab() + xlim(0,2)
```
