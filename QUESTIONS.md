# QUESTIONS

## DCA  
- 论文地址: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6344535/pdf/41467_2018_Article_7931.pdf
- 论文补充材料：https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-07931-2/MediaObjects/41467_2018_7931_MOESM1_ESM.pdf
- 论文数据集: 见论文Data Availability Statement (对应数据集已经给出了网址)
- GITHUB源码：https://github.com/theislab/dca
- 任务1:实现DCA算法。
- 任务2:如论文中所写，在模拟数据集上，如论文所写分两类的模拟数据和六类的模拟数据这两种情况测试您自行实现的DCA的效果。
- 任务3:如论文中所写，在Zheng数据集上测试DCA是否能保留聚类的结构，并且测试对于论文中几个marker gene(某一类细胞特有的表达的基因)您自行实现的DCA是否能够恢复其在特定聚类中的表达。
- 任务4:如论文中所写，测试您自行实现的DCA是否能将连续分化的血细胞嵌入一条按照分化顺序决定的曲线上。
- 任务5:如论文中所写，测试您自行实现的DCA是否能在组织RNA测序和单细胞RNA测序两个数据集的对照下恢复单细胞RNA测序的一些基因表达。
- 任务6:如论文中所写，测试您自行实现的DCA是否能恢复转录组表达和细胞表面标志物蛋白之间的相关性。

## SC3  
- 论文地址:https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5410170/
- 论文nature地址: https://www.nature.com/articles/nmeth.4236#Ack1
- 论文数据集:https://hemberg-lab.github.io/scRNA.seq.datasets/ 在此地址下搜索对应数据集名称即可
- GITHUB源码：[SC3](https://github.com/hemberg-lab/SC3) and [FastSC3](https://github.com/hemberg-lab/FastSC3) and [SC3s](https://github.com/hemberg-lab/sc3s) and [论文图表生成源码](https://github.com/hemberg-lab/SC3-paper-figures)
- [R package](http://bioconductor.org/packages/release/bioc/html/SC3.html) and [SC3 R tutorial](http://bioconductor.org/packages/release/bioc/vignettes/SC3/inst/doc/SC3.html)
- 任务1: 自行实现SC3算法(不要求SVM优化)，并在论文中给出的12个数据集上进行测试。
- 任务2: 在R语言环境下安装SC3的软件包,将问题1中的效果和该软件本身的效果进行比较。
- 任务3: 对于你实现的SC3算法进行SVM优化，并在论文中提到的大规模数据集或者在其它大规模数据集(大于5000个，比如上述地址中的Baron)中进行验证，并和SC3使用SVM优化的效果进行比较。
- 任务4: 注意该论文中的表述存在一些错误，其中提到对于距离矩阵使用PCA，然而其实际实现时并为将每个样本对应的向量投影到PCA的前k个主成成分所对应的向量上，而是直接使用前k个主成成分作为PCA的结果。这有一点像PCoA(主坐标分析),然而对于PCoA来说，结果并未在每个维度乘上对应的特征值。或许这种处理也和KPCA(核技巧在PCA上的应用)有所联系。写一写你对上述题干的理解，并以理论或实际有进一步表现为目标尝试给出一个改进。更进一步，你是否能尝试改变SC3的各个距离度量模块以及降维模块，以获得更好的聚类效果。

