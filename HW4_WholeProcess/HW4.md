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

## 5. Постройте ML-дерево в RAxML-NG, используя выбранную модель. (Приведите код или описание действий.)

```
raxml-ng --check --msa SUP35_aln_prank.trim.fas  --model TIM3+G4 --prefix SUP35_raxml_test
raxml-ng --msa SUP35_aln_prank.trim.fas --model TIM3+G4 --prefix SUP35_raxml --threads 2 --seed 222 
```
## 6. Отрисуйте полученное дерево (лучшее ML-дерево) и покажите рисунок.
+1 балл за скрипт на любом языке, который принимает на вход название файла и автоматически рисует качественное дерево (скрипт можно привести внутри отчёта, вряд ли он очень длинный).

## 7. Как выбрать модель в ModelFinder (можно через IQ-TREE)?

## 8. Какая модель эволюции была признана наиболее подходящей для нашего выравнивания?
(Скопируйте нужную часть вывода программы).
Отличаются ли модели, выбранные ModelTest и ModelFinder, и насколько сильно?

## 9. Постройте ML-дерево в IQ-TREE, используя выбранную модель. (Приведите код или описание действий.)

## 10. Отрисуйте полученное дерево (лучшее ML-дерево) и покажите рисунок.