# coding4bioinfo
coding test for bioinfo, bioinfomatics.

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

### 关于ZINB分布
- [负二项分布](https://zh.wikipedia.org/wiki/%E8%B4%9F%E4%BA%8C%E9%A1%B9%E5%88%86%E5%B8%83)
- [负二项分布及其应用](https://zhuanlan.zhihu.com/p/111632687)
- [**零膨胀负二项模型**](https://www.jianshu.com/p/149ff509fe7f): 从测序的机理上建模scRNA-seq数据的分布，进而假设泊松分布，该分布均值方差相等，但是实际上有over-dispersion问题，因此尝试对泊松分布的`lambda`加gamma分布的先验，两个汇总推导得到X服从负二项分布。此时一个问题是，technical noises and biological variability (or dropout)导致零值非常多，这是经典的zero-inflation问题。于是进一步得到ZINB分布模型，其中新增的π可以视为真实的基因表达值被观测为0的概率。同时，ZINB分布也可以reparameteried为基于mean和dispersion parameter的构造。
- [基于gamma分布先验的泊松分布Poisson–Gamma Mixture等价于负二项分布推导证明](https://gregorygundersen.com/blog/2019/09/16/poisson-gamma-nb/)
- [单细胞RNA-seq数据分布的选择](https://zhuanlan.zhihu.com/p/95299303)
----------------------

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

### PCA降维再理解
- [PCA的数学原理](http://blog.codinglabs.org/articles/pca-tutorial.html)
- [延申阅读KPCA](https://blog.csdn.net/lyn5284767/article/details/81509059)
- [参考晓雯的整理资料，关于距离度量和降维方法](https://github.com/sxwenny/job/blob/master/%E6%9C%BA%E5%99%A8%E5%AD%A6%E4%B9%A0.md)
- 一些后续可以尝试的思路。线性降维：SVD、非线性降维：KPCA、LLE(Local linearly embedding)、autoencoder and tSNE。

----------------------


## Daily Schedule
- 3.15 通读DCA文章，理解所解决的问题：scRNAseq的denoise和imputation。通过autoencoder，既可以做到dimension reduction又可以做到imputation；深度理解DCA的算法原理，思路比较清晰；暂对ZINB分布模型有一定疑惑。
- 3.16 继续理解ZINB分布，以及相应的NB分布，zero-inflation模型，NB分布的两种起源形式（一种基于概率分布的定义，另一种基于泊松分布+gamma分布的推导）和ZINB的两种参数化形式。
- 3.16 通读SC3文章，理解所解决的问题：scRNAseq的cell cluster问题。SC3算法更像是一套cluster workflow not algorithm，总结而言：特征选择（基因过滤）、距离矩阵计算（三种距离形式）、矩阵变换（PCA和laplacian）、kmeans for multi combination、consensus cluster and hierarchical cluster. workflow又有一些扩展：当样本量过大时，即cell num>5000，此时先采样5000样本用workflow进行unsupervised cluster，然后得到伪标签用SVM进行supervised learning提高speed. 目前原理理解下来，整个思路比较清晰。
- 3.17 回顾PCA原理，考虑task4的可能的改进方案：关于距离度量和降维。是否可以在降维这里用上autoencoder？
- 3.17 [安装R环境](https://blog.csdn.net/Joshua_HIT/article/details/73741139) [VScode配置R](https://blog.csdn.net/u011262253/article/details/113837720) 
- 3.18 学习R的基本知识 [SC3分析博客](http://t.zoukankan.com/leezx-p-10878506.html) 尝试基于SC3demo跑通实验流程，多次遇到因为内存和CPU挂掉的问题，目前该SC3.R代码已经可以正常走通。但仅在20000*90的scRNAseq数据上。大数据量上该咋办？耗时也较长。
- 3.19 
  - 数据的normalization[参考](http://www.360doc.com/content/18/0112/02/50153987_721216719.shtml)：对read counts采用基因长度[对于某一个基因而言]和测序深度[对于某一个样本而言]
  - 下载数据 
  - [python的scanpy库](https://scanpy.readthedocs.io/en/latest/index.html)以及pip直接安装，[报错 ImportError: DLL load failed while importing utilsextension；解决方案：单独重装tables==3.6.1即可](https://github.com/theislab/scanpy/issues/2108)
  - python sc3 流程大框架撰写
  - [距离度量](https://cloud.tencent.com/developer/article/1406436)