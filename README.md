# coding4bioinfo
coding test for bioinfo, bioinfomatics.

| ARI  (est_k)          | Yan        | Biase    | Goolam   |  Deng1 (FIXED) |  Deng2 (FIXED) | Klein      | Baron      |
| ------------          | ---------- |----------|----------|----------      | ----------     | ---------- | ---------- |
| Rcode                 | 0.6584(6)  | 0.9871(5)| 0.5973(6)| ~~0.4268(8)~~ 0.4439(9)| 0.7876(9)| 0.6178(12)| 0.2977(37)|
| Python(base)          | 0.6584(6)  | 0.9871(5)| 0.5973(6)| ~~0.3625(9)~~ 0.3827(9)| 0.6648(9)| 0.7291(12)| 0.3276(37)|
|+ normdists            | 0.6584(6)  | 0.9871(5)| 0.5973(6)|               0.3647(9)| 0.6557(9)| 0.6870(12)| -         |
|- laplacian            | 0.7212(6)  | 0.9839(5)| 0.5973(6)|               0.4333(9)| 0.7745(9)| 0.8224(12)| 0.4179(37)|
|- laplacian pca->decom | 0.6584(6)  | 0.9839(5)| 0.5973(6)|               0.4642(9)| 0.8176(9)| 0.5312(12)| - |
|- laplacian + normdists| 0.6584(6)  | 0.9394(5)| 0.5860(6)|               0.4504(9)| 0.7975(9)| 0.6141(12)| 0.2983(37)|
|- laplacian + JS dist  | 0.7080(6)  | 0.9839(5)| 0.5973(6)|               0.3023(9)| 0.5858(9)| 0.8070(12)| - |
|- laplacian + UMAP*    | 0.7902(6)  | 0.8302(5)| 0.5927(6)|               0.3606(9)| 0.6616(9)| 0.5548(12)| - |
|- laplacian - spearman | 0.7212(6)  | 0.9839(5)| 0.5973(6)|               0.4440(9)| 0.7888(9)| 0.8247(12)| 0.3757(37)|
|- laplacian - spearman +rbf|0.7212(6)|0.9871(5)|0.5927(6)|0.4721(9)|0.8366(9)|0.8167(12)|0.4228(37)|
|-laplacian -spearman +linear|0.7212(6) |0.9807(5)|0.5973(6)|0.4343(9)|0.7722(9)|0.8372(12)|0.4208(37)|
|-laplacian -spearman +rbf pca>decom (similarity)|- |-|-|-|0.6870(9)|0.3402(12)|-|


## Update
- 3.30 完成[阶段性报告输出](http://htmlpreview.github.io/?https://github.com/gatsby2016/coding4bioinfo/blob/main/docs/SC3-exp-report.html)；更新部分数据集cluster图和rbf的h5ad文件；Baron和klein文件太大未上传.

- 3.29 针对PCoA和kPCA的等价推导，有以下疑问：
  - 能否对输出的dist matrix直接应用PCoA？
  - 能否直接对dist matrix做特征分解？此时相当于把dist matrix当成kernel matrix
  - 能否修改dist matrix，即，改变不同的kernel function得到不同的kernel matrix？有方案：RBFkernel 或JS散度 或Jaccard系数度量
  - 接第一条PCoA（MDS），思考降维能否采用embedding的方式实现？即，IsoMap, LLE, spectral embedding, tSNE or UMAP？在以上输入dissimilarity matrix的情况下，仅spectral embedding, tSNE or UMAP可用。
  - LLE延伸：laplacian eigenmap，取最小的d个非零特征值对应的特征向量；进而延申为，谱聚类spectral cluster方法。 继续延申，距离度量，是否可以采用k个最近点的方式？
- 3.28
  - 主要完成workflow与PCoA及KPCA的关联[推导证明](docs/proof.pdf)
  - 基于推导结果验证修改距离度量为RBF kernel计算，并直接进行特征分解，结果较差，待更新
- 3.27
  - 完成Rcode中Baron数据svm mixed model实验验证对比；
  - 验证dist calculation中remove spearman dist的结果；【分析：spearman相关系数对分量排序离散化，可以理解为筛选当前cell的marker gene 同时确定pairwise的cell是否也具有同样趋势；】
  - 距离度量部分：是否需要修改现有euclidean距离度量为标准化euclidean距离，以消除不同维度分量量纲影响；【分析：修改现有euclidean为标准化没有意义，因为对于某个样本来说，某个gene可能的确是marker gene，进而导致在该样本下gene expression value higher but lower in other cells. 同时其他gene的expression lower。】
  - 把kmeans结果之后CSPA之前的multi-res直接ARI分析。【**发现**：laplacian采用后k个特征值对应的特征向量作为降维进行kmeans结果ARI和pca降维kmeans结果相似；而采用前k个特征值对应的特征向量的kmeans后几乎都预测为一类，导致ARI很低！所以，尝试将laplacian降维关闭！结果见表！】
  - 距离度量添加JS散度距离；结果有下降，耗时巨大；
  - 降维部分增加基于优化的非线性降维UMAP算法；优化学习由于收敛问题引入了不确定性；
  - 距离度量后，Pearson和spearman距离范围都是[0,2],但是欧式距离没有上限；已进行归一化。欧式距离先对样本单位化再归一化。结果没有提升；
  - 推导原方案与kernelPCA的关联；
    - [proof1](docs/proof1.png)
    - [proof2](docs/proof2.png)
- 3.26 
  - 新增支持Klein数据；
  - python脚本优化计算spearman和pearson距离效率；
  - 将Rcode结果送入到k_cluster_matrix中进行后续操作，结果和Rcode完全一致；证实整个python workflow仅在kmeans聚类有差异；
  - run **baron**数据集，目前整个runtime控制在1h以内（57min）；排查结果差的原因：estimate k cluster (37) > true cluster (14)，手动指定k值为14 ARI为0.56
  - 重构SC3脚本，剥离SC3类文件与main文件；增加对数据变量adata结构的h5ad格式保存与载入，保存位置为`SC3/Results`
- 3.25 
  - 修复Baron数据上获取`n_dim`时的dtype错误；
  - 修复gene_filter操作与R不一致问题；原有方式基于count进行filter，然后应用在logcount上，已对齐；
  - 继续验证Rcode与Python不一致问题；已确定是kmeans算法问题；Rcode使用`Hartigan-Wong`算法[ref](https://doi.org/10.2307/2346830) or [here](https://sci-hub.se/https://doi.org/10.2307/2346830)；Python使用`Elkan algorithm with random(or kmeans++)`；两者样本迭代过程有差异。[排查确定方式](SC3/Python/SC3.py#L303)
- 3.24 
  - 新增支持Deng数据集；已完成结果比对，~~但结果和文章不一致，**待排查**。~~ 
  - **排查与文章不一致**结论：[数据网站](https://hemberg-lab.github.io/scRNA.seq.datasets/mouse/edev/#deng)上写的deng数据是`22431 features and 268 samples`；但是下载下来按照其bash脚本处理后，`22958 features and 267 samples`。
  - **排查Rcode和python不一致**结论：Rcode经过了一次duplicated后特征维度降为7861维。原因：缺失`feature_symbol`信息。
  - **问题已修复！！！** [数据网站给的deng.sh生成gene_name.txt时1这里少加了域号$导致txt为空没有维度信息进而导致一系列错误。](SC3/Data/Deng/deng.sh#L53) 现已和`22431 features and 268 samples`数据对齐。
  - **修复Goolam数据上的稍微差异问题 [On CPM calculation](SC3/Python/SC3.py#L61)**
  - [调整kmeans聚类](SC3/Python/SC3.py#L27) `Number of time the k-means algorithm will be run with different centroid seeds` 
  - 排查确定**kmeans随机初始化影响**较大；已对齐kmeans输入的情况下，经过kmeans聚类结果不管是python内部还是python与Rcode都有较大不同；目前暂无法解决该问题。**具体地：deng数据上** 
    - 在python内部，kmeans的`init`参数指定`random`与默认会导致ARI 0.6648与0.5553的差别；其中，指定`random`主要为了对齐Rcode；
    - python与Rcode之间，同样`random centroids from row of data`的情况下，会存在0.5553与0.7876的差异。
- 3.23 新增Goolam数据集及data process R脚本，修复create_sce脚本calculateCPM传参问题；同步更新SC3的R和python脚本支持Goolam数据；目前已完成结果比对，见表。
  - understanding PCoA [PCoA主坐标分析;](https://qinqianshan.com/math/gradient/pcoa-analysis/) [principal_coordinates_analysis;](https://www.sequentix.de/gelquest/help/principal_coordinates_analysis.htm) [group discussion;](https://groups.google.com/g/qiime-forum/c/i-2uhMk-Lug)  [PCoA and NMDS;](https://www.davidzeleny.net/anadat-r/doku.php/en:pcoa_nmds)
  - [聚类算法中不同相似度相异度矩阵计算](https://blog.csdn.net/enochzhu/article/details/109769648)
  - 重新理解SC3设计流程。第一步特征选择filter;第二步计算distance matrix；此时用到了euclidean and pearson and spearman；三种度量，一种距离，两种相似性系数，而euclidean和pearson都是基于连续量计算，spearman则对变量进行rank，也就变成了顺序的、离散量计算；第二步完成后得到的其实是dissimilarity matrix；然后第三步，对dissimilarity matrix应用PCA和graph laplacian进行降维，说是PCA，但其实又不是PCA，第一输入不是data matrix而是dissimilarity matrix，第二不是用前K个主成分变换，而是取得前K个主成分；而取前K个最大特征值对应特征向量的操作和PCoA很像，而且PCoA也应用在dissimilarity matrix。**PCoA analysis is also called metric multi-dimensional scaling**。第四步，对降维（或者说变换）后的数据进行多重kmeans聚类，并进行Similarity matrix融合，最后进行层次聚类获得final result。
- 3.22 1. 新增支持R code based SC3 workflow定量结果输出，直接run脚本即可。实现Rcode和python的横向比对。2. 增加层次聚类后一致性矩阵图绘制。
- 3.21 新增支持基于SVM-mixed hybrid model 并且重构SC3类
  - hybrid model `50 cluster+40 svm`：`Yan:{'6': 0.6913063345361797}` 耗时也有近半下降
  - hybrid model `40 cluster+ remained svm`： `{'5': 0.9871491842026757}`效果未改善
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
- 3.19 
  - 数据的normalization[参考](http://www.360doc.com/content/18/0112/02/50153987_721216719.shtml)：对read counts采用基因长度[对于某一个基因而言]和测序深度[对于某一个样本而言]
  - 下载数据 
  - [python的scanpy库](https://scanpy.readthedocs.io/en/latest/index.html)以及pip直接安装，[报错 ImportError: DLL load failed while importing utilsextension；解决方案：tables==3.6.1单独重装](https://github.com/theislab/scanpy/issues/2108)
  - python sc3 流程大框架撰写
  - [距离度量](https://cloud.tencent.com/developer/article/1406436)- 3.19 
  - 数据的normalization[参考](http://www.360doc.com/content/18/0112/02/50153987_721216719.shtml)：对read counts采用基因长度[对于某一个基因而言]和测序深度[对于某一个样本而言]
  - 下载数据 
  - [python的scanpy库](https://scanpy.readthedocs.io/en/latest/index.html)以及pip直接安装，[报错 ImportError: DLL load failed while importing utilsextension；解决方案：tables==3.6.1单独重装](https://github.com/theislab/scanpy/issues/2108)
  - python sc3 流程大框架撰写
  - [距离度量](https://cloud.tencent.com/developer/article/1406436)
- 3.18 学习R的基本知识 [SC3分析博客](http://t.zoukankan.com/leezx-p-10878506.html) 尝试基于SC3demo跑通实验流程，多次遇到因为内存和CPU挂掉的问题，目前该SC3.R代码已经可以正常走通。但仅在20000*90的scRNAseq数据上。大数据量上该咋办？耗时也较长。
- 3.17 [安装R环境](https://blog.csdn.net/Joshua_HIT/article/details/73741139) [VScode配置R](https://blog.csdn.net/u011262253/article/details/113837720) 
- 3.17 回顾PCA原理，考虑task4的可能的改进方案：关于距离度量和降维。是否可以在降维这里用上autoencoder？
- 3.16 通读SC3文章，理解所解决的问题：scRNAseq的cell cluster问题。SC3算法更像是一套cluster workflow not algorithm，总结而言：特征选择（基因过滤）、距离矩阵计算（三种距离形式）、矩阵变换（PCA和laplacian）、kmeans for multi combination、consensus cluster and hierarchical cluster. workflow又有一些扩展：当样本量过大时，即cell num>5000，此时先采样5000样本用workflow进行unsupervised cluster，然后得到伪标签用SVM进行supervised learning提高speed. 目前原理理解下来，整个思路比较清晰。
- 3.16 继续理解ZINB分布，以及相应的NB分布，zero-inflation模型，NB分布的两种起源形式（一种基于概率分布的定义，另一种基于泊松分布+gamma分布的推导）和ZINB的两种参数化形式。
- 3.15 通读DCA文章，理解所解决的问题：scRNAseq的denoise和imputation。通过autoencoder，既可以做到dimension reduction又可以做到imputation；深度理解DCA的算法原理，思路比较清晰；暂对ZINB分布模型有一定疑惑。


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
- 一些后续可以尝试的思路。特征分解角度；学习embedding（非线性）角度：KPCA、LLE(Local linearly embedding)、autoencoder and tSNE， UMAP。
- [Ten quick tips for effective dimensionality reduction](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006907) and [相异性矩阵降维](https://baijiahao.baidu.com/s?id=1638747337747160251&wfr=spider&for=pc)
- [Paper:Evaluation of UMAP as an alternative to t-SNE for single-cell data](https://www.biorxiv.org/content/10.1101/298430v1.full)
- [表征图数据，图embedding](https://baijiahao.baidu.com/s?id=1672090485285983725&wfr=spider&for=pc)
- [拉普拉斯谱聚类](https://cloud.tencent.com/developer/article/1155839)
- [从PCA到PCoA到KPCA](https://blog.csdn.net/janehong1314/article/details/86618311) and [PCoA分析](https://www.cnblogs.com/lyon2014/p/3991058.html)
- [矩阵获取距离度量矩阵](https://blog.csdn.net/frankzd/article/details/80251042)
- [LLE原理及推导过程](https://blog.csdn.net/scott198510/article/details/76099630)
----------------------

## INSTALL (to be updated)
- [python conda env](sc3_env.yaml)
----------------

## TBD
- ~~python hybrid model with SVM training support~~
- ~~R code ARI quantitative metric completely compared to Python's~~
- ~~python plots func support for visualization, such as `sc3_plot_consensus`~~ and `sc3_plot_expression`
- ~~spearman distance slow calculation [尝试自己实现python版本和c++版，随着维度上升，效率十分慢于pandas；已优化，替换pandas计算为scipy.stats.spearmanr]~~
- ~~more datasets support [add deng, goolam, klein and baron dataset, others cannot downloaded from hemberg-lab site]~~
- ~~dist mearsurement and dims reduction support~~
- further analysis implementation on the cluster results

-------------

## 工程优化考虑
- 特征选择 [特征选择策略]
- 距离度量 []
- 特征分解进行降维 []
- kmeans聚类 []
- CSPA combination [assign weights]
- 层次聚类 