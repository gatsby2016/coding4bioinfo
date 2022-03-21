# coding4bioinfo
coding test for bioinfo, bioinfomatics.


## Update
- 3.21 新增支持基于SVM-mixed hybrid model 并且重构SC3类
  - hybrid model `50 cluster+40 svm`：`Yan:{'6': 0.6913063345361797}` 耗时也有近半下降
  - hybrid model `40 cluster+ remained svm`： `{'5': 0.9871491842026757}`效果未改善

-----------------

## DCA 关于ZINB分布
- [负二项分布](https://zh.wikipedia.org/wiki/%E8%B4%9F%E4%BA%8C%E9%A1%B9%E5%88%86%E5%B8%83)
- [负二项分布及其应用](https://zhuanlan.zhihu.com/p/111632687)
- [**零膨胀负二项模型**](https://www.jianshu.com/p/149ff509fe7f): 从测序的机理上建模scRNA-seq数据的分布，进而假设泊松分布，该分布均值方差相等，但是实际上有over-dispersion问题，因此尝试对泊松分布的`lambda`加gamma分布的先验，两个汇总推导得到X服从负二项分布。此时一个问题是，technical noises and biological variability (or dropout)导致零值非常多，这是经典的zero-inflation问题。于是进一步得到ZINB分布模型，其中新增的π可以视为真实的基因表达值被观测为0的概率。同时，ZINB分布也可以reparameteried为基于mean和dispersion parameter的构造。
- [基于gamma分布先验的泊松分布Poisson–Gamma Mixture等价于负二项分布推导证明](https://gregorygundersen.com/blog/2019/09/16/poisson-gamma-nb/)
- [单细胞RNA-seq数据分布的选择](https://zhuanlan.zhihu.com/p/95299303)
----------------------

## SC3 距离度量降维与聚类
- [PCA的数学原理](http://blog.codinglabs.org/articles/pca-tutorial.html)
- [延申阅读KPCA](https://blog.csdn.net/lyn5284767/article/details/81509059)
- [理解不同的距离度量和降维方法](https://github.com/sxwenny/job/blob/master/%E6%9C%BA%E5%99%A8%E5%AD%A6%E4%B9%A0.md)
- 一些后续可以尝试的思路。线性降维：SVD、非线性降维：KPCA、LLE(Local linearly embedding)、autoencoder and tSNE。

----------------------


## Daily Work
- 3.15 通读DCA文章，理解所解决的问题：scRNAseq的denoise和imputation。通过autoencoder，既可以做到dimension reduction又可以做到imputation；深度理解DCA的算法原理，思路比较清晰；暂对ZINB分布模型有一定疑惑。
- 3.16 继续理解ZINB分布，以及相应的NB分布，zero-inflation模型，NB分布的两种起源形式（一种基于概率分布的定义，另一种基于泊松分布+gamma分布的推导）和ZINB的两种参数化形式。
- 3.16 通读SC3文章，理解所解决的问题：scRNAseq的cell cluster问题。SC3算法更像是一套cluster workflow not algorithm，总结而言：特征选择（基因过滤）、距离矩阵计算（三种距离形式）、矩阵变换（PCA和laplacian）、kmeans for multi combination、consensus cluster and hierarchical cluster. workflow又有一些扩展：当样本量过大时，即cell num>5000，此时先采样5000样本用workflow进行unsupervised cluster，然后得到伪标签用SVM进行supervised learning提高speed. 目前原理理解下来，整个思路比较清晰。
- 3.17 回顾PCA原理，考虑task4的可能的改进方案：关于距离度量和降维。是否可以在降维这里用上autoencoder？
- 3.17 [安装R环境](https://blog.csdn.net/Joshua_HIT/article/details/73741139) [VScode配置R](https://blog.csdn.net/u011262253/article/details/113837720) 
- 3.18 学习R的基本知识 [SC3分析博客](http://t.zoukankan.com/leezx-p-10878506.html) 尝试基于SC3demo跑通实验流程，多次遇到因为内存和CPU挂掉的问题，目前该SC3.R代码已经可以正常走通。但仅在20000*90的scRNAseq数据上。大数据量上该咋办？耗时也较长。
- 3.19 
  - 数据的normalization[参考](http://www.360doc.com/content/18/0112/02/50153987_721216719.shtml)：对read counts采用基因长度[对于某一个基因而言]和测序深度[对于某一个样本而言]
  - 下载数据 
  - [python的scanpy库](https://scanpy.readthedocs.io/en/latest/index.html)以及pip直接安装，[报错 ImportError: DLL load failed while importing utilsextension；解决方案：tables==3.6.1单独重装](https://github.com/theislab/scanpy/issues/2108)
  - python sc3 流程大框架撰写
  - [距离度量](https://cloud.tencent.com/developer/article/1406436)
- 3.20 完成**SC3 basic framework** (only cluster without SVM training) python版本开发。数据处理，三个距离矩阵计算，两个降维算法实现并且都已经和Rcode一致对齐。后续的kmeans聚类，转换相似矩阵，层次聚类等方法输出结果未和Rcode比较，仅观察最终输出ARI定量结果(下方)，目前应该没有问题。该版本已经在biase和yan数据上跑通。
  - Yan数据 `>>> adjustedRandIndex {'6': 0.6584306303041297}`, 
  - biase数据 `>>> adjustedRandIndex {'5': 0.9871491842026757}`
  - [理解距离度量：Pearson euclidean and cosine](https://blog.csdn.net/sixtyfour/article/details/80354164)
  - [机器学习常用距离定义与计算1](https://zhuanlan.zhihu.com/p/101277851) [2](https://zhuanlan.zhihu.com/p/266490448)
  - [Pearson VS Spearman](https://blog.csdn.net/lambsnow/article/details/79972145)
  - Pearson： $1 - \frac{(u - \bar{u}) \cdot (v - \bar{v})}
                  {{||(u - \bar{u})||}_2 {||(v - \bar{v})||}_2}$
  - [decomposition.PCA](https://www.cnblogs.com/pinard/p/6243025.html)
  - [graph laplacian matrix](https://zhuanlan.zhihu.com/p/25096844)
  - [VScode R debug](https://blog.csdn.net/qq_42679415/article/details/120374896)
  - [因为graph laplacian没对齐 尝试调试看R调用C函数的norm_laplacian](https://www.cnblogs.com/lotusto/p/5740297.html#:~:text=R%E8%AF%AD%E8%A8%80%E8%B0%83%E7%94%A8C%2B%2B%20Rcpp%E5%8C%85%E6%98%AF%E4%B8%80%E4%B8%AA%E6%89%93%E9%80%9AR%E8%AF%AD%E8%A8%80%E5%92%8CC%2B%2B%E8%AF%AD%E8%A8%80%E7%9A%84%E9%80%9A%E4%BF%A1%E7%BB%84%E4%BB%B6%E5%8C%85%EF%BC%8C%E6%8F%90%E4%BE%9B%E4%BA%86R%E8%AF%AD%E8%A8%80%E5%92%8CC%2B%2B%E5%87%BD%E6%95%B0%E7%9A%84%E7%9B%B8%E4%BA%92%E8%B0%83%E7%94%A8%E3%80%82,R%E8%AF%AD%E8%A8%80%E5%92%8CC%2B%2B%E8%AF%AD%E8%A8%80%E7%9A%84%E6%95%B0%E6%8D%AE%E7%B1%BB%E5%9E%8B%E9%80%9A%E8%BF%87Rcpp%E5%8C%85%E8%BF%9B%E8%A1%8C%E5%AE%8C%E6%95%B4%E7%9A%84%E6%98%A0%E5%B0%84%E3%80%82%20R%E8%AF%AD%E8%A8%80%E8%B7%A8%E7%95%8C%E8%B0%83%E7%94%A8C%2B%2B) [rtools4.0官方](https://cran.r-project.org/bin/windows/Rtools/rtools40.html) [利用Rcpp和RcppArmadillo创建R包](https://blog.csdn.net/iamsuperman2/article/details/77103568) 
  - [搞懂norm_laplacian代码 基于官方RcppArmadillo的c++ API](http://arma.sourceforge.net/docs.html#each_colrow)
- 3.21 新增支持基于SVM-mixed hybrid model 并且重构SC3类
  - hybrid model `50 cluster+40 svm`：`Yan:{'6': 0.6913063345361797}` 耗时也有近半下降
  - hybrid model `40 cluster+ remained svm`： `{'5': 0.9871491842026757}`

----------------
## INSTALL (to be updated)
- [python conda env](sc3_env.yaml)
----------------

## TBD
- ~~python hybrid model with SVM training support~~
- python plots func support for visualization, such as `sc3_plot_consensus` and `sc3_plot_expression`
- spearman distance slow calculation
- more datasets support
- dist mearsurement and dims reduction support
