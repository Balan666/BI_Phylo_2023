# Python: Bio::Phylo

```
import requests
import matplotlib
import random
from io import StringIO
from Bio import Phylo
```

##### 11. читает дерево https://www.jasondavies.com/tree-of-life/life.txt 
`raw_tree = StringIO(requests.get('https://www.jasondavies.com/tree-of-life/life.txt').text)`

##### 12. рисует это дерево псевдографикой (draw_ascii);
```
tree1 = Phylo.read(raw_tree, "newick")
Phylo.draw_ascii(tree1)
```

##### 13. рисует дерево с помощью draw (картинка должна быть внутри скомпилированного документа);
`Phylo.draw(tree1)`

##### 14. сохраняет изображение дерева в растровый формат (png) и векторный (svg/pdf);
```
matplotlib.pyplot.savefig("py_tree1_phylo.svg",  branch_labels=lambda x: x.branch_length)
matplotlib.pyplot.savefig("py_tree1_phylo.png")
```

##### 15. рисует дерево в читаемом виде (размер шрифта)
```
#Phylo.draw_graphviz(tree1) #networkx

tree1.clade[0, 1].color = "blue"
matplotlib.rc('font', size=1)
matplotlib.pyplot.figure(figsize=(24,12))
Phylo.draw(tree1, do_show = False)
matplotlib.pyplot.savefig("py_tree1_phylo_blue.png", dpi=600)
```

## Python. ETE (ETE3)

##### 16. читает дерево;
`raw_tree = requests.get('https://www.jasondavies.com/tree-of-life/life.txt').text`

##### 17. рисует это дерево;
`tree2 = Tree(raw_tree, format=1)`

##### 18. рисует это дерево в читаемом виде;
`tree2.render("py_tree2_ete3.pdf")`

В формате круга:

```
circular_style = TreeStyle()
circular_style.mode = "c"
circular_style.scale = 20
tree2.render("py_tree2_ete3_circ.pdf", tree_style=circular_style)
```

##### 19. вырезает (функция prune) из дерева случайный набор из 42 листьев;

```
ss = random.sample(tree2.get_leaf_names(), 42)
tree2.prune(ss)
```

##### 20. рисует обрезанное дерево.

`tree2.render("py_tree2_ete3_random.pdf")`

##### TreeStyle from ETE3
```
ts = TreeStyle()
ts.branch_vertical_margin = 15

nstyle = NodeStyle()
nstyle["shape"] = "square"
nstyle["size"] = 5 
nstyle["fgcolor"] = "lightblue"
for leaf in tree2.traverse(): 
   leaf.set_style(nstyle)

tree2.render("py_tree2_ete3_random_improved.png", tree_style=ts)
```
