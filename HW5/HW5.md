## 1. версии использованных программ.
1. IQ-Tree >=2.1.3 (http://www.iqtree.org/, если не взлетает, ищите внизу странички web service или https://anaconda.org/bioconda/iqtree).
Обратите внимание на версию. Текущая = 2.2.0, проще всего поставить её.

пакеты для R: phytools, ape, ggtree, ggplot, ggpubr.

```{r}
library(easypackages)
install_packages(c('phytools', 'ape', 'ggtree', 'ggplot2', 'ggpubr'))
libraries(c('phytools', 'ape', 'ggtree', 'ggplot2', 'ggpubr'))
```

## 2. BOOTSTREP

Базовая команда:
```
iqtree2 -s SUP35_aln_prank.trim.fas -m TIM3+F+G4 -pre SUP35_TIM3
```
>Оцениваем устойчивость топологии:
>как запустить построение дерева в iqtree2, но
>с генерацией 100 реплик обычного бутстрепа?

чтобы было быстро, поставим 10 повторов
```
time iqtree2 -s SUP35_aln_prank.trim.fas -m TIM3+F+G4 -redo -pre SUP35_TIM3_b -b 10
real	0m10,024s
```
но вообще надо минимум 100
```
time iqtree2 -s SUP35_aln_prank.trim.fas -m TIM3+F+G4 -redo -pre SUP35_TIM3_b -b 100
```
## 3 ULTRAFAST BOOTSTREP
ультрабыстрый бутстреп, как мне показалось, раз в 10 быстрее

ультрабыстрый бутстреп дает значения % больше, его нужно бырать строже, по 95%
```
time iqtree2 -s SUP35_aln_prank.trim.fas -m TIM3+F+G4 -redo -pre SUP35_TIM3_ufb -bb 1000
```
## 4. approximate bayes test
```
iqtree2 -s SUP35_aln_prank.trim.fas -m TIM3+F+G4 -pre SUP35_TIM3_B_alrt_abayes -bb 1000 -alrt 1000 -abayes
0m1,744s
```
#https://itol.embl.de/ 

Поддержка: */bias/bootstrep

```{r}
library(ggtree)
tree_alrt_abayes_ufb <- read.tree("SUP35_TIM3_B_alrt_abayes.treefile")
ggtree(tree_alrt_abayes_ufb) + 
  geom_tiplab() + geom_nodelab() +
  geom_treescale() + xlim(0, 0.7)
# funny labels
label <- tree_alrt_abayes_ufb$node.label
alrt <- as.numeric(sapply(strsplit(label, "/"), FUN = "[", 1)) #sun
abayes <- as.numeric(sapply(strsplit(label, "/"), FUN = "[", 2)) #yin yang
ufb <- as.numeric(sapply(strsplit(label, "/"), FUN = "[", 3)) #star
large_alrt <- ifelse(alrt > 70, intToUtf8(9728), "")
large_abayes <- ifelse(abayes > 0.7, intToUtf8(9775), "")
large_ufb <- ifelse(ufb > 95, intToUtf8(9733), "")
newlabel <- paste0(large_alrt, large_abayes, large_ufb)
tree_alrt_abayes_ufb$node.label <- newlabel
ggtree(tree_alrt_abayes_ufb) + 
  geom_tiplab() + geom_nodelab(nudge_x = -.01, nudge_y = .1) +
  geom_treescale() + xlim(0, 0.7)
``` 

<img src='https://github.com/Balan666/BI_Phylo_2023/blob/main/HW5/1a723d2a-d801-4397-8329-fdaf9789c2e7.png?raw=true'>
<img src='https://github.com/Balan666/BI_Phylo_2023/blob/main/HW5/88739747-41c4-4b33-8eda-7b7acc64ea98.png?raw=true'>

## Корни
можно брать за корень:
## 5.
###### внешнюю группу (sufficiently close), на таксономический уровень выше
```
iqtree2 -s SUP35_aln_prank.trim.fas -m TIM3+F+G4 -pre SUP35_TIM3_root_outgroup -bb 1000 -alrt 1000 -abayes  -o SUP35_Kla_AB039749,SUP35_Agos_ATCC_10895_NM_211584
```
(аргумент -о)

ITOL:
пикча
## 6.
###### midpoint rooting

#https://github.com/mooreryan/midpoint-root
```{r}
#install.packages("phytools")
library(phytools)
midpoint.root(tree_alrt_abayes_ufb)
```

ITOL:
<img src='https://github.com/Balan666/BI_Phylo_2023/blob/main/HW5/midpointroot.png?raw=true'>

## 7.
###### Необратимые модели
учитыват течение времени

```
iqtree2 -s SUP35_aln_prank.trim.fas -m TIM3+F+G4 -pre SUP35_TIM3_root_auto --model-joint 12.12 -B 1000
```

выдает nexus дерево и rootsrap (точность оставляет желать лучшего)

with a certain model:
```
iqtree2 -s SUP35_aln_prank.trim.fas -m JC -pre SUP35_JC -bb 1000 -alrt 1000 -abayes -o SUP35_Kla_AB039749,SUP35_Agos_ATCC_10895_NM_211584
```
## 8. Дерево с поддержкой корня (rootstrap)
<img src='https://github.com/Balan666/BI_Phylo_2023/blob/main/HW5/rootstraptree.png?raw=true'>

## 9.То же самое дерево с использованием заведомо неправильной модели JC+G4 и поддержками
<img src='https://github.com/Balan666/BI_Phylo_2023/blob/main/HW5/wrong.png?raw=true'>

## 10. Деревья, построенные с моделями TIM3+F+G4 и JC (SUP35_JC.treefile и SUP35_TIM3_root_outgroup.treefile).
выдает два дерева. сравнивать их можно много где

#http://phylo.io/ 

#https://beta.phylo.io/viewer/#


```{r}
library(ggtree)
treeTIM3 <- read.tree("SUP35_TIM3_root_outgroup.treefile")
treeJC <- read.tree("SUP35_JC.treefile")
library(ggplot2)
tim3tree <- 
  ggtree(treeTIM3) + geom_tiplab() +
  geom_nodelab(color = "blue", nudge_x = -.05) + 
  theme_tree2() + 
  xlim(0,1) + ggtitle("TIM3")
jctree <- 
  ggtree(treeJC) + geom_tiplab() +
  geom_nodelab(color = "red", nudge_x = -.05) + 
  theme_tree2() + 
  xlim(0,1) + ggtitle("JC")
library(ggpubr)
ggarrange(tim3tree, jctree)
association <- cbind(treeTIM3$tip.label, treeJC$tip.label)
cophyloplot(treeTIM3, treeJC, assoc=association, length.line=4, space=28, gap=3)
library(phytools)
trees.cophylo<-cophylo(treeTIM3, treeJC, rotate = TRUE)
png("cophylo.png", width = 1200, height = 800)
plot(trees.cophylo, link.type="curved",link.lwd=4,
     link.lty="solid",link.col=make.transparent("red", 0.25), size = 1)
dev.off()
```
ITOL:
<img src='https://github.com/Balan666/BI_Phylo_2023/blob/main/HW5/compare.png?raw=true'>
