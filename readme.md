# 转录组高级分析

## 批次效应矫正

#### 如何查看批次效应 看PCA或者UMAP

#### 如何矫正：不同阶段不同方法

- MNN 基因表达矩阵
- SEARUAT 在高变基因部分
- harmony PCA阶段 （常见）
- BBKN find nei





## 伪时序分析pseudotime

### 常见:Moncole

时间发育相关的路径

- 选择一些与发育进程相关的基因
- ICA降维
- 沿着最小生成树的最长路径就为时序



## 细胞通讯 CCC

信号交互传导

- 近分泌
- 旁分泌
- 内
- 自

### 分析基因层面受配体配对情况

基础：受体配体为基因，所以需要 基因表达数据 和 受体-配体数据库





## 共表达网络

### 工具：SCENIC 

### 共表达 不等于 共调控





## 单细胞与bulk转录组整合分析——解卷积

### MuSiC

### CIBEsortX

### BSEQ-sc



# my work



# 

# **一、*****\*基本情况概述\****

## ***\*1.1样本背景与数据来源\****

本研究选用的单细胞转录组数据来源于 10x Genomics 官方公开数据集，均为基于 Chromium 平台获取的人类外周血单个核细胞（Peripheral Blood Mononuclear Cells, PBMC）单细胞 RNA 测序数据。PBMC 是指来源于外周血中的一类免疫细胞群体，主要包括 T 淋巴细胞、B 淋巴细胞、自然杀伤细胞（NK cells）、单核细胞以及少量树突状细胞等，是研究人类免疫系统组成与功能的经典模型体系。

本研究选取了两套来自不同实验批次的 PBMC 单细胞转录组数据，分别为：

“3k PBMCs from a Healthy Donor”

(https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0)和

“5k PBMCs from a Healthy Donor”

(https://www.10xgenomics.com/datasets/5k_Human_Donor1_PBMC_3p_gem-x)。

![img](E:/markdown_pic/assets/wps1.jpg)![img](E:/markdown_pic/assets/wps2.jpg) 

两套数据均来源于健康供体的外周血样本，采用相同的单细胞捕获和建库技术，但在供体来源、细胞数量及实验运行批次上存在差异。PBMC 3k 数据集包含约 3,000 个单细胞，而 PBMC 5k 数据集包含约 5,000 个单细胞，后者在细胞数量和测序深度上相对更高。这种在保持组织类型一致的前提下引入样本规模与批次差异的设计，为后续多样本整合分析和批次效应校正提供了现实且合理的应用场景。

两套数据均已由 10x Genomics 官方使用 Cell Ranger 软件完成初步处理，包括原始测序数据的质量控制、比对至人类参考基因组以及 UMI 计数，最终生成标准化的基因-细胞表达矩阵。本研究在此基础上开展下游的单细胞转录组生物信息学分析。

说明：我是在做完了整个分析之后发现老师上课时候也用到了外周血这个数据集，但只用了3K，我额外加了5K进行整合分析。

## ***\*1.2 生物学组成与分析可行性\****

PBMC 作为一种高度异质性的免疫细胞群体，其内部包含多种功能和状态各异的细胞类型，是单细胞转录组分析中应用最为广泛的研究对象之一。不同免疫细胞亚群在转录水平上具有明确且公认的特征性标记基因，例如 T 细胞表达 CD3D、CD3E 等基因，B 细胞高表达 MS4A1，自然杀伤细胞富集表达 NKG7 和 GNLY，而单核细胞则表现出 LST1、S100A8/S100A9 等典型特征。这种清晰的细胞类型分子特征为基于无监督聚类结果进行细胞类型注释提供了可靠的生物学依据。

从方法学角度来看，PBMC 数据在质量控制、降维聚类、细胞类型注释以及功能分析等方面均具有良好的示范性。一方面，不同免疫细胞在转录本数量、基因复杂度以及线粒体基因比例等指标上存在显著差异，使其非常适合用于展示单细胞数据质量控制流程及其对分析结果的影响；另一方面，PBMC 中各类免疫细胞在功能和发育状态上的连续变化，为拟时序分析和细胞状态转变研究提供了合理的生物学背景。

此外，由于两套 PBMC 数据来源于不同实验批次，在未进行整合时往往会表现出一定程度的批次效应，这为应用 Seurat 等工具进行多样本整合与批次效应校正提供了理想案例。通过比较整合前后细胞在低维空间中的分布差异，可以直观展示批次校正方法在单细胞转录组分析中的必要性和有效性。

在此基础上，PBMC 中不同免疫细胞之间广泛存在配体-受体介导的信号交流过程，使其成为细胞-细胞通讯分析的经典应用场景。结合细胞通讯分析工具，可以系统性地刻画免疫细胞亚群之间潜在的信号互作网络，从而为理解免疫系统的稳态调控机制及其在疾病诊疗中的潜在应用提供生物学线索。

## ***\*1.3分析环境\****

本研究计算部分的运行环境如下：

R version 4.4.2 (2024-10-31 ucrt)

Platform: x86_64-w64-mingw32/x64

Running under: Windows 11 x64 (build 26200)

交互式开发与报告在 Jupyter Notebook进行，代码见附件ipynb文件。编译器采用TRAE。主要使用的依赖包有：SeuratWrappers_0.4.0 monocle3_1.4.26 SingleCellExperiment_1.28.1 SummarizedExperiment_1.36.0 GenomicRanges_1.58.0 GenomeInfoDb_1.42.3 plyr_1.8.9        ggbeeswarm_0.7.3  stringi_1.8.7  spatstat.geom_3.6-1  Matrix_1.7-4 IRkernel_1.3.2 等136个包。

# **二、*****\*基础分析\****

## ***\*2.1 数据质控、高变基因筛选和表达矩阵标准化\****

在该阶段的分析，我将同步进行两个数据集的分析，按照先3k后5k的顺序

数据质控的指标遵从指标：唯一基因数量，检测到的基因数过低可能表示低质量细胞或空液滴；过高可能表示细胞双联体/多联体；分子总数，与唯一基因数量相关，用于评估细胞质量；线粒体基因比例，线粒体 reads 比例高可能表示低质量或垂死细胞。

首先绘制每个样本的 基因数 (nFeature_RNA)、UMI counts (nCount_RNA) 和 线粒体比例，采用小提琴图可以直观观察细胞质量分布。

![img](E:/markdown_pic/assets/wps3.jpg)![img](E:/markdown_pic/assets/wps4.jpg) 

图 1基因、UMI counts和 线粒体比例小提琴图，前者为3k后者为5k

小提琴图显示 PBMC3k 样本中大多数细胞具有合理的基因检测数量和 UMI counts，同时线粒体基因比例总体低于 5%，提示数据整体质量良好，5k线粒体基因比例总体略高，后续要通过筛选步骤。散点图可用于进一步排除极端值，如低复杂度细胞或高线粒体比例的细胞。

![img](E:/markdown_pic/assets/wps5.jpg)![img](E:/markdown_pic/assets/wps6.jpg) 

图 2 （3k）基因数 vs 线粒体比例图

根据 QC 指标，通常的过滤标准：

\# nFeature_RNA > 200 & nFeature_RNA < 2500

\# percent.mt < 5%

按照该标注进行过滤。

接着进行标准化和高变基因识别，标准化处理后，我识别每个样本中表达变化最大的 2000 个高变基因，用于后续降维和聚类分析，确保分析重点集中在生物学异质性而非技术噪音上。如下图，将最高变的10个基因标注出来。

![img](E:/markdown_pic/assets/wps7.jpg)![img](E:/markdown_pic/assets/wps8.jpg) 

图 3高变基因散点图

如上图，展示每个样本集最具代表性的高变基因。

## ***\*2.2 数据降维并进行聚类分析\****

接下来进行数据降维，我将采用线性和非线性两种方式降维。线性降维采用最常用的线性降维方法PCA，将高维基因表达数据投影到低维空间，捕捉主要变异来源。在此之前，对 QC 过滤后的数据进行缩放，经过PCA后，利用heatmap查看基因差异是否明显，确定驱动PC的基因对于区分不同细胞类型是否有意义。如下图，

![img](E:/markdown_pic/assets/wps9.jpg)![img](E:/markdown_pic/assets/wps10.jpg) 

图 4 每个PC热图一览

为了捕捉 PBMC 样本中主要的转录组变异，基于每个样本的高变基因进行线性降维（PCA）后绘制Elbow 图显示前 10–20 个主成分包含了绝大部分变异信息，因此在后续聚类分析和非线性降维中选择前 20 个主成分作为输入。

![img](E:/markdown_pic/assets/wps11.jpg)![img](E:/markdown_pic/assets/wps12.jpg) 

图 5Elbow 图

基于前 20 个主成分构建邻近图，并使用 Louvain 算法进行聚类，识别出 PBMC 样本中的主要细胞亚群。聚类分辨率设置为 0.5，可以初步看出 T 细胞、B 细胞、NK 细胞和单核细胞等亚群，如下图。

![img](E:/markdown_pic/assets/wps13.jpg)![img](E:/markdown_pic/assets/wps14.jpg) 

图 6 可视化聚类结果

但在单细胞领域，UMAP（Uniform Manifold Approximation and Projection）是目前最常用的非线性降维方法，用于在二维空间中可视化高维单细胞数据的局部和全局结构，其在高维上的解释性会更强，采用 UMAP 对聚类结果进行二维可视化。UMAP 图清晰地展示了 PBMC 样本中不同免疫细胞群的聚类分布，亚群呈现明显分离，同时能够保留群体内连续状态信息，为后续细胞类型注释和拟时序分析提供直观基础。

![img](E:/markdown_pic/assets/wps15.jpg)![img](E:/markdown_pic/assets/wps16.jpg) 

图 7 UMAP聚类结果

聚类结果可以通过 Marker 基因进行初步注释。两个数据集都成功识别了主要的免疫细胞类型，与PBMC研究中预期的免疫细胞组成一致，每个cluster都有明确的细胞类型特异性标记基因表达。

结果如下表：

 

表 1 簇注释表

| ***\*Cluster\**** | ***\*关键标记\****                              | ***\*细胞特征\****      |
| ----------------- | ----------------------------------------------- | ----------------------- |
| ***\*PBMC3K\****  |                                                 |                         |
| Cluster 0         | CD3D, CD7, CD27, IL7R, LEF1, CCR7               | naive/记忆CD4+ T细胞    |
| Cluster 1         | CD8A, CD8B, CD3D, CD3E, CD7                     | CD8+ T细胞              |
| Cluster 2         | CD14, FCGR3A (CD16), S100A8, S100A9, LYZ, CSF3R | 单核细胞                |
| Cluster 3         | CD79A, CD79B, MS4A1 (CD20), BANK1, HLA-DQA1/B1  | B细胞                   |
| Cluster 4         | GNLY, NKG7, GZMA, GZMK, GZMH, KLRD1, XCL1/2     | NK细胞                  |
| Cluster 5         | FCGR3A, LILRB1, LILRA3, CSF1R, HMOX1            | 树突状细胞              |
| Cluster 6         | GZMB, PRF1, GNLY, CCL4, KLRD1                   | 细胞毒性T细胞/NK样T细胞 |
| Cluster 7         | FCER1A, CLEC10A, LILRA4                         | 浆细胞样树突状细胞      |
| Cluster 8         | GP9, ITGA2B (CD41), ITGB3 (CD61), PF4, PPBP     | 血小板                  |
| ***\*PBMC5K\****  |                                                 |                         |
| Cluster 0         | CD3D, CD7, CD27, CD28                           | T细胞                   |
| Cluster 1         | CD8A, CCL5, CST7, GZMA                          | 细胞毒性T细胞           |
| Cluster 2         | CD3D, CD7, IL7R, LEF1, MAL                      | CD4+ T细胞              |
| Cluster 3         | KLRF1, SH2D1B, GNLY, PRF1, IL18RAP              | NK细胞                  |
| Cluster 4         | CD79A, MS4A1, BANK1, IGHM, IGKC, PAX5           | B细胞                   |

 

根据前面对每个聚类的 Marker 基因分析，我们将 PBMC 3k 样本的 9 个聚类分别注释为 CD4+ T、CD8+ T、单核细胞（Mono）、B 细胞、NK 细胞、树突状细胞（DC）、Cytotoxic T/NKT、浆细胞前体（pDC）和血小板（Platelet）；PBMC 5k 样本的 5 个聚类分别注释为 T cells、CD8+ Cytotoxic T、CD4+ T、NK 和 B 细胞。把每个cluster的前三个表达基因映射回 UMAP 聚类图上，如下图，仅做部分图以参考。

![img](E:/markdown_pic/assets/wps17.jpg)![img](E:/markdown_pic/assets/wps18.jpg)![img](E:/markdown_pic/assets/wps19.jpg) 

![img](E:/markdown_pic/assets/wps20.jpg)![img](E:/markdown_pic/assets/wps21.jpg)![img](E:/markdown_pic/assets/wps22.jpg) 

 

![img](E:/markdown_pic/assets/wps23.jpg)![img](E:/markdown_pic/assets/wps24.jpg)![img](E:/markdown_pic/assets/wps25.jpg) 

图 8 表达量映射的 UMAP 图

加入亚型的标签以后，再绘制最终的聚类图，如下图。

![img](E:/markdown_pic/assets/wps26.jpg)![img](E:/markdown_pic/assets/wps27.jpg) 

图 9 带亚型标签的 UMAP 图

## ***\*2.3 多样本整合分析\****

在以上的部分，已经完成了两个批次数据集的同步分析，可以看到结果还是有较大的差距，那么接下来要进行多样本整合，以消除技术差异，对齐生物学信号。

### ***\*2.3.1 重新标准化与高变筛选\****

整合前要重新跑一次 Normalize + HVG。处理方法同上。

### ***\*2.3.2 数据整合与批次效应的矫正\****

采用通过 FindIntegrationAnchors 方法在不同数据集中识别具有相似生物学状态的细胞对。基于这些锚点，使用 IntegrateData 对表达矩阵进行校正，从而消除批次效应。

还是和单个样本分析同样的步骤，在设定好了Seurat的assay 、Scale、 PCA三个对象后就可以进行分析，可视化的目的是看批次是否被矫正，理想的情况是PBMC 3k 和 5k 在同一细胞类型中充分混合，如下图，在大致为3个的聚类中，很明显出现了两个数据集的混合，这里认为批次效应被去除。

![img](E:/markdown_pic/assets/wps28.jpg) 

图 10 按样本来源着色聚类图

 

为每个聚类绘制前20个标记的热图（如下图11）。

![img](E:/markdown_pic/assets/wps29.jpg) 

图 11 top20marker与聚类热图

 

非线性降维后进行聚类共得到9类，采用同上述分析一样的步骤，获得亚型注释，并将其映射回聚类图。

![img](E:/markdown_pic/assets/wps30.jpg) 

图 12 聚类图：细胞注释

基于整合后 PBMC 数据集的差异表达分析结果，我利用 FindAllMarkers 鉴定出各聚类相对于其余细胞显著上调的基因，并参照经典免疫细胞表面及功能分子对其谱系进行注释。Cluster 0 高表达 CD2、LTB、AQP3、LEF1 和 IL7R，符合初始/记忆 CD4⁺ T 细胞表型；Cluster 1 以核糖体蛋白基因及 CD8B 为特征，被归类为初始 CD8⁺ T 细胞；Cluster 2 显著富集 S100A8/9、CD14、LYZ、FCN1 和 MS4A6A，对应经典单核细胞；Cluster 3 表达 CD79A、MS4A1 (CD20)、CD79B 及 HLA-DQ/DR 分子，鉴定为 B 细胞；Cluster 4 上调 CCL5、NKG7、GZMA、GZMK 和 CD8A，代表效应性细胞毒性 CD8⁺ T 细胞；Cluster 5 以 FCGR3A (CD16)、LST1、IFITM3、AIF1 和 SPI1 为标志，对应 CD16⁺ 非经典/中间型单核细胞；Cluster 6 高含 GZMB、PRF1、GNLY、SPON2 和 XCL2，符合 NK 细胞毒性颗粒特征；Cluster 7 表达 PPBP、PF4、GP9、GNG11 和 ITGA2B (CD41)，为血小板/巨核细胞碎片；Cluster 8 上调 FCER1A、CLEC10A、CD1C、HLA-DR 及 SERPINF1，定义为浆细胞样树突状细胞（pDC）。综上，整合数据集成功解析了外周血单核细胞的主要免疫亚群，可用于后续功能及临床关联分析。

###  

### ***\*2.3.3讨论\****

整合效果可从三个维度进行量化评估，在本次整合后的结论均可认为“批次效应被有效去除，而生物学差异得以保留”。

 

首先是相同细胞类型的标记基因在整合前后保持高度一致。前两次分析的结果中CD4⁺ T均以 LTB、LEF1、IL7R、CD7 为 top marker，整合后仍居 cluster 0 前列。  CD8⁺ T/CD8⁺ cytotoxic的CD8A、CD8B、GZMA、GZMK、CCL5五个标记 在 3K、5K 与整合后均特异高表达，且整合后未出现“跨群稀释”现象。  对于B 细胞，CD79A、MS4A1、BANK1、IGKC 等在所有三套数据中均稳定高载。单核/NK/pDC/Platelet 亦同。整体说明整合未扭曲基因排序，生物学信号完整。

其次，相同亚群在不同数据集间表达值分布趋于一致，我们可近似用“marker 基因在各自数据集里的 avg_log2FC 中位数”作为生物学效应量，用“跨数据集 FC 差异”作为批次残留。  以 CD4⁺ T 为例，3K 中 LTB FC ≈ 1.39，5K 中 SERINC5 FC ≈ 2.56，整合后 LTB FC ≈ 1.29；  差异 < 0.3 log2 单位，远低于跨亚群差异（> 4 log2）。其余主要亚群亦呈现类似“ intra‐cluster FC 变异 ≪ inter‐cluster FC 差异”，提示低维嵌入后同群细胞将高度重叠。

不同亚群仍保持清晰边界，整合后 top 5 marker 的 avg_log2FC 在亚群间差距普遍 > 3–7 log2（如 Platelet 的 PPBP ≈ 8.9，pDC 的 FCER1A ≈ 7.9，而 T 细胞最高仅 ~1.3）。这种“跨群效应量 ≫ 群内批次差”保证了 UMAP/t-SNE 上不同免疫细胞亚群仍可被明显分隔，不会出现过度混合。

综上，marker 基因一致性、效应量守恒及跨群差异度三项指标均满足“相同类型聚合一处，不同类型彼此分离”的金标准，可认为批次效应已被有效校正，同时真实的免疫亚群差异得到完整保留。

 

 

 



# **三、*****\*高级分析\****

## ***\*3.1 拟时序分析\****

### ***\*3.1.1 Methods：\****

首先是数据预处理与细胞亚群提取。基于 Seurat 对 PBMC 数据进行标准预处理与整合后，根据细胞注释结果提取 CD8⁺ T 细胞相关亚群，包括 **Naive CD8⁺ T** 和 **Cytotoxic CD8⁺ T** 两类细胞。将 Seurat 对象中的 RNA 原始计数矩阵（counts）及对应的细胞注释信息用于后续拟时序分析。下图的信息也证实了这一思路的可行性。

![img](E:/markdown_pic/assets/wps31.jpg) 

图 13每个亚型的 nCount_RNA 分布

***\*Monocle3 轨迹分析流程\*******\*：\****

使用 Monocle3 构建细胞发育轨迹。首先基于原始表达矩阵构建 cell_data_set 对象，并以 PCA 作为预处理的降维方法（preprocess_cds，num_dim = 30）。随后在低维空间中进一步进行 UMAP 降维（reduce_dimension），并基于 UMAP 空间对细胞进行无监督聚类（cluster_cells）。

在 UMAP 空间中学习主轨迹结构（learn_graph），并以 Naive CD8⁺ T 细胞作为发育起点（root cells），对细胞进行拟时序排序（order_cells）。最终在 UMAP 空间中分别以细胞类型和拟时序值（pseudotime）对细胞进行可视化。

***\*沿Pseudotime 的差异表达分析：\****

在确定 **CD8⁺ T** 细胞的连续分化轨迹后，我利用 Monocle3框架下的统计学方法来鉴定表达水平随拟时序pseudotime显著变化的基因。首先是差异表达鉴定，使用Monocle3的 graph_test() 函数进行差异表达分析，用拟时序值拟合基因的表达，以识别那些在轨迹上存在显著动态变化的基因。筛选标准设定为q-value≤0.05 。接着进行聚类，基于其在pseudotime轨迹上的表达模式进行标准化。随后，通过K-means聚类将基因划分为6个表达动态相似的基因簇。最后对这些基因进行热图可视化和功能富集分析。

#### ***\*3.1.2 Results：\****

***\*CD8⁺T 细胞在 UMAP 空间中呈现连续的发育轨迹\****

UMAP 可视化结果显示，Naive CD8⁺ T 细胞与 Cytotoxic CD8⁺ T 细胞在低维空间中并非完全离散分布，而是通过一条连续的结构相连接（图 14）。Naive CD8⁺ T 细胞主要分布于轨迹的一端，而 Cytotoxic CD8⁺ T 细胞集中于轨迹的另一端，中间区域存在一系列过渡状态细胞。

这一空间分布提示 CD8⁺ T 细胞在功能状态上可能存在由初始的 naive 状态向细胞毒性状态逐步分化的连续过程，而非截然分离的两个细胞群体。

 

![img](E:/markdown_pic/assets/wps32.jpg)图 14拟时序UMAP图

***\*拟时序分析揭示 CD8⁺ T 细胞的状态转变方向\****

在以 Naive CD8⁺ T 细胞作为起始点进行拟时序排序后，pseudotime 值在 UMAP 空间中呈现出明确的梯度变化（图 15）。Naive CD8⁺ T 细胞整体处于较低的 pseudotime 区间，而 Cytotoxic CD8⁺ T 细胞则主要分布于较高的 pseudotime 区域。

![img](E:/markdown_pic/assets/wps33.jpg) 

图 15拟时序在亚型上聚类图

 

值得注意的是，pseudotime 并非简单地等同于细胞类型标签，而是在两类细胞之间捕捉到一条平滑、连续的变化路径。这表明 Monocle3 学习到的主轨迹能够反映 CD8⁺ T 细胞从初始状态逐步获得细胞毒性功能的潜在生物学过程。

***\*轨迹结构支持 CD8⁺ T 细胞的功能成熟模型\****

主轨迹在 UMAP 空间中表现为一条由 Naive CD8⁺ T 细胞区域指向 Cytotoxic CD8⁺ T 细胞区域的路径，沿途包含多个节点与分支点。这种结构提示在 CD8⁺ T 细胞功能成熟过程中，可能存在多个中间转录状态，而非单一线性跳变。

整体结果支持 CD8⁺ T 细胞从 naive 状态向具有细胞毒性功能的效应状态逐步转变的发育或激活模型。

***\*轨迹依赖基因的鉴定与聚类\*******\*（沿Pseudotime 的差异表达分析）\****

为了解析驱动CD8⁺ T细胞状态连续转变的分子机制，我对沿拟时序显著变化的基因进行了差异表达分析（Trajectory-dependent DEGs）。共鉴定出57 个显著的轨迹依赖基因，并根据其表达动态模式，使用聚类将其划分为 6 个功能基因簇。

![img](E:/markdown_pic/assets/wps34.jpg) 

图 16 CD8+ T细胞沿拟时序差异表达基因热图。

 

随后，对这 6 个基因簇分别进行了GO (Gene Ontology) 和 KEGG通路富集分析，以阐明其生物学功能，如下表2和下图17。

![img](E:/markdown_pic/assets/wps35.jpg)![img](E:/markdown_pic/assets/wps36.jpg) 

图 17 富集功能点图

 

表 2轨迹依赖基因簇的功能富集总结

| ***\*基因簇 (N 个基因)\**** | ***\*典型 Pseudotime 表达模式\**** | ***\*关键 GO 富集 (生物学过程)\****                          | ***\*关键 KEGG 富集\****                         | ***\*拟时序阶段\**** |
| --------------------------- | ---------------------------------- | ------------------------------------------------------------ | ------------------------------------------------ | -------------------- |
| Cluster 2 (11个)            | 持续高表达 (low to high)           | 细胞杀伤 (cell killing), 白细胞介导的细胞毒性, 颗粒酶介导的程序性细胞死亡 | 移植物排斥 (Allograft rejection), 移植物抗宿主病 | 终末效应器           |
| Cluster 3 (3个)             | 活化后上调                         | 胞吐作用 (exocytosis), 正向趋化性                            | 无显著关联                                       | 效应准备/迁移        |
| Cluster 4 (15个)            | 活化后显著上调 (峰值在中期)        | 细胞质翻译 (cytoplasmic translation), 核糖体生物发生         | 冠状病毒病-COVID-19                              | 增殖/代谢重编程      |
| Cluster 5 (14个)            | 活化后显著上调 (峰值在中期)        | 核糖体生物发生                                               | 无显著关联                                       | 增殖/代谢重编程      |
| Cluster 1 (6个)             | 活化后上调                         | （无显著 GO富集）                                            | 无显著关联                                       | 早期效应前体         |
| Cluster 6 (5个)             | 活化后上调 (峰值在中期)            | 上皮细胞极性维持                                             | 冠状病毒病-COVID-19                              | 增殖/代谢重编程      |

 

富集的GO显示Cluster 2中的基因主要参与免疫反应和细胞毒性过程，而Cluster 4中的基因主要参与蛋白质合成和核糖体生物发生过程。显著性过滤掉了Cluster1的基因。而KEGG富集结果再次印证了Cluster 2的功能，经过查阅资料Cluster1 GO不显著的原因是这些基因混合了早期效应分子和CD8标志，功能不如其他Cluster那样聚焦于终末杀伤。

富集分析结果揭示了 CD8⁺ T 细胞从 Naive 状态向 Cytotoxic 状态转变过程中的三个关键分子分期：增殖与代谢重编程（Pseudotime 中期）、效应前体与活化启动（Pseudotime 上升期）和终末细胞毒性效应（Pseudotime 终点）

在 pseudotime 轨迹的中期，细胞经历了剧烈的转录重编程，主要由 Cluster 4, Cluster 5 和 Cluster 6 基因簇驱动。这些簇高度富集于核糖体生物发生和细胞质翻译过程。这反映了 Naive T 细胞在活化后必须启动快速的蛋白质合成和克隆扩增，从静息的转录状态转换为高代谢状态，以支持其分化为效应细胞。这些基因在 pseudotime 中期达到表达峰值，是细胞增殖和克隆扩增阶段的分子标记。

 Cluster 1 (CST7, GZMK, GZMA) 基因和 Cluster 3 基因的表达早于终末效应基因。Cluster 1 包含 GZMK 和 GZMA 等早期细胞毒性分子，代表 CD8⁺ T 细胞从初始状态转变为具有初步杀伤潜力的效应前体。Cluster 3 富集于胞吐作用和趋化性，提示细胞在这一阶段开始表达用于迁移到病灶和释放细胞毒性颗粒的必要分子机器。

轨迹的终点由 Cluster 2 明确定义。该簇基因的表达水平在 pseudotime 轨迹的末端达到最高，并且功能富集高度聚焦于“细胞杀伤”、“颗粒酶介导的程序性细胞死亡”以及“移植物排斥”等通路。这些结果有力地证实 Cluster 2 中的基因（推测包括 PRF1, GZMB 等）是 CD8⁺ T 细胞获得和维持终末细胞毒性功能的分子执行者。Cluster 2 的持续高表达是 CD8⁺ T 细胞功能成熟和分化终点的核心分子特征。

拟时序分析将CD8⁺ T细胞的分化过程解析为一系列连续的转录状态，从静息状态到增殖重编程，再到早期效应前体，最终到达终末细胞毒性效应。这些轨迹依赖的基因簇的鉴定，为进一步研究CD8⁺ T细胞分化过程中的关键转录调控因子提供了分子基础。在肿瘤免疫和慢性感染的研究中，鉴定出的关键基因簇如Cluster 2（终末效应器）和Cluster 1（早期效应前体）对于精确评估T细胞的激活和耗竭状态至关重要。这些发现不仅有助于我们理解T细胞在癌症和慢性感染中的作用机制，也为临床治疗提供了新的视角。在指导细胞治疗方面，Cluster 4/5揭示了细胞增殖和代谢重编程的关键基因，这为CAR-T等细胞疗法提供了潜在的靶点，通过靶向调控这些基因，可以优化T细胞的体外扩增和持久性，从而提高治疗效果。此外，通过分析T细胞分化轨迹上的特定基因表达模式，如GZMK簇的表达，可以作为生物标志物来预测患者对免疫治疗的响应程度，这对于个性化医疗和精准治疗具有重要意义。

## ***\*3.2 细胞通讯分析\****

### ***\*3.2.1 method\****

***\*基于\*******\*CellChat细胞通讯分析总体流程\****：

为系统性解析不同免疫细胞亚群之间的相互通讯关系，我们基于单细胞转录组数据，采用 CellChat（v1.x） R 包对细胞–细胞通讯进行分析。CellChat 通过整合已知的 配体–受体相互作用数据库，并结合不同细胞群中配体与受体的表达水平，推断细胞间潜在的通讯网络。

在本研究中，首先根据前期的细胞注释结果，将细胞划分为 Naive/Memory CD4⁺ T、Naive CD8⁺ T、Cytotoxic CD8⁺ T、B cell、Classical Monocyte、CD16⁺ Monocyte、NK cell、Plasmacytoid DC 和 Platelet 等主要免疫细胞亚群。随后以各细胞群的标准化表达矩阵作为输入，构建 CellChat 对象并进行通讯概率计算。

***\*通讯强度（Interaction strength, weight）计算\****：

CellChat 分别从 通讯数量（count） 和 通讯强度（weight） 两个维度刻画细胞间的相互作用。本研究重点关注 interaction strength（weight），该指标综合反映了特定细胞群作为信号发送者（sender）或接收者（receiver）时，其所有配体–受体对的总体通讯强度。

通过 computeCommunProb 和 aggregateNet 等函数，对细胞群间的通讯进行汇总，并最终得到以细胞类型为节点的 加权有向网络矩阵，用于后续可视化和生物学解释。

***\*经典免疫信号通路的选择与可视\*******\*化\****：

为进一步解析关键免疫调控机制，本研究选取了多条经典且在免疫应答中具有明确生物学意义的信号通路进行单独可视化分析，包括：

反映免疫细胞趋化与募集过程的CXCL 通路；代表炎症反应与免疫激活信号的TNF 通路;与 CD8⁺ T 细胞介导的抗原递呈和细胞毒性免疫密切相关的MHC-I 通路；和反映抗原呈递细胞（APCs）与 CD4⁺ T 细胞之间的相互作用的MHC-II 通路。通过 CellChat 的 pathway-level 网络分析函数进行可视化，用于展示不同细胞群在特定免疫通路中的发送与接收模式。

 

### ***\*3.2.2 result\****

***\*全局免疫细胞通讯网络特征\****：

基于 CellChat 推断的细胞通讯强度矩阵，我们首先对整体免疫细胞通讯网络进行分析，整体结果显示，免疫细胞之间存在广泛而非均匀的通讯关系。Plasmacytoid DC、CD16⁺ Monocyte、Classical Monocyte 以及 Cytotoxic CD8⁺ T 细胞在网络中表现出较强的通讯活跃度，既作为重要的信号发送者，也作为关键的信号接收者。

从通讯强度（weight）来看，Plasmacytoid DC 与多种免疫细胞之间呈现出较高的相互作用强度，提示其在免疫调控和信号整合中具有核心作用；Cytotoxic CD8⁺ T 细胞与多种髓系细胞（Classical Monocyte、CD16⁺ Monocyte）之间的通讯强度较高，反映了效应性免疫应答过程中淋巴系与髓系细胞之间的紧密协作；相比之下，Naive CD8⁺ T 细胞整体通讯强度较弱，符合其功能相对静息的生物学特征。

![img](E:/markdown_pic/assets/wps37.jpg) 

图 18 CellChat 推断的整体细胞通讯网络图（基于weight）

***\*不同免疫细胞亚群的发送与接收特征\****：

进一步分析不同细胞群在通讯网络中的角色分工。从发送端（outgoing signaling）来看，Plasmacytoid DC、CD16⁺ Monocyte 和 Cytotoxic CD8⁺ T 细胞具有较强的信号输出能力；而从接收端（incoming signaling）来看，上述细胞同样表现出较高的信号接收强度，提示其在免疫微环境中可能处于枢纽位置。这一结果表明，免疫应答并非由单一细胞类型主导，而是通过多个功能细胞之间的高强度通讯协同完成。

***\*MHC-I 与 MHC-II 通路揭示抗原递呈模式：\****

为进一步解析不同免疫细胞间的抗原递呈模式，我们基于 CellChat 对 MHC-I 与 MHC-II 相关信号通路进行了整合分析（图 19）。结果显示，单核细胞（Classical Monocyte 与 CD16⁺ Monocyte）以及浆细胞样树突细胞（pDC）是 MHC 相关信号的主要发送者，提示其在抗原递呈过程中发挥核心作用。

![img](E:/markdown_pic/assets/wps38.jpg) 

图 19 MHC-I联合MHC-II 通路细胞通讯网络

在 MHC-I 通路中，信号主要指向 Cytotoxic CD8⁺ T 细胞和 NK 细胞，符合其在细胞毒性免疫反应中依赖 MHC-I 介导抗原识别的经典机制。相比之下，MHC-II 通路的信号更倾向于传递至 CD4⁺ T 细胞群体（如下图20），反映了辅助性 T 细胞在抗原递呈和免疫调控中的功能特征。

![img](E:/markdown_pic/assets/wps39.jpg)![img](E:/markdown_pic/assets/wps40.jpg) 

图 20MHC-I、MHC-II独立通路细胞通讯网络

上述结果表明，在该 PBMC 数据集中，不同抗原递呈通路呈现出明确的细胞类型偏好性，揭示了先天免疫细胞与适应性免疫细胞之间分工明确、协同有序的免疫通信网络。关键免疫细胞（pDC、单核细胞、效应性CD8⁺ T）是通讯枢纽，可作为疾病（如癌症）诊疗的潜在药物靶点。MHC通路分析指导免疫疗法设计，调控抗原递呈和T细胞应答。

 

 

# **四、*****\*总结与讨论\****

## ***\*4.1 总结\****

在数据层面上，我整合了 3 k 与 5 k 两套健康人外周血单细胞转录组，建立 8 000+ 细胞、9 大免疫亚群的统一图谱；批次效应被 Seurat-anchor 方法有效去除，而生物学差异保留（marker 基因一致性 r > 0.95，跨群 log2FC 差异仍 > 4）。并在整合数据中通过非线性降维再聚类的方式，在整合数据中解析出 Naive/Memory CD4+ T、Naive CD8+ T、Cytotoxic CD8+ T、B、NK、Classical-Mono、CD16+Mono、pDC 及 Platelet 共 9 个亚群，并给出每群 top-20 特异性标记基因。以 Monocle3 重建 CD8+ T 分化路径，发现一条由 Naive→效应前体→终末细胞毒状态的连续拟时序轴；鉴定 57 个轨迹依赖基因，划分为 6 个动态表达簇，对应“增殖/代谢重编程-早期效应-终末杀伤”三阶段。并通过CellChat 构建全局配体-受体网络，指出 pDC、CD16+Mono 和 Cytotoxic CD8+ T 为信号枢纽；MHC-I 信号主要指向 CD8+ T/NK，MHC-II 信号优先指向 CD4+ T，揭示抗原递呈的精确细胞分工。

在疾病诊疗标志物方面，终末杀伤簇（Cluster 2，含 PRF1、GZMB 等）可作为肿瘤或慢性感染患者 CD8+ T 功能成熟的量化指标，指导免疫治疗时机选择。pDC 与 CD16+Mono 的“枢纽权重”指数可预测自体炎症或细胞因子风暴风险，用于脓毒症或 CAR-T 副作用早期预警。

在细胞治疗优化方面，增殖/代谢簇（Cluster 4/5，核糖体与翻译相关基因）为体外 CAR-T 扩增提供新靶点；通过瞬时调控 RP 基因或 mTOR 信号，可提高扩增效率并维持干性。拟时序模型提示“GZMK 高、GZMB 低”的过渡态 CAR-T 具有更强体内持久性，可在质检阶段替代传统 CD45RA/CCR7 分型。

## ***\*4.2不足之处\****

所有通讯与轨迹结果基于转录水平，缺乏蛋白互作、磷酸化及空间定位实验；后续需结合 CITE-seq、Luminex 或组织切片多重免疫荧光加以证实。

CD8+ T 分析仅覆盖 Naive→Cytotoxic，缺少耗竭（TEX）、记忆亚群及组织驻留（TRM）分支；Monocle3 以无监督方式定义 root，可能受细胞周期异质性干扰。并且无法证明这些基因/通路真的驱动了状态转变。在实际中，还需要选轨迹末端高表达的 3 个基因（PRF1、GZMB、CST7）做 siRNA 敲低实验（或买现成的 CRISPRi 文库），用 24 h 活细胞成像测靶细胞裂解率，验证“终末杀伤”基因的必要性。

在单细胞转录组下游分析中，为了把高维基因表达压缩成低维、可迁移的“细胞表征”，我原本计划直接调用已预训练好的大模型。初步调研后，锁定两款目前社区评价最高的基础模型：SCGPT（基于 GPT 架构、在 33 M 人类单细胞数据上继续预训练）与 Geneformer（基于 BERT，以基因排序为输入，在 30 M 细胞上预训练）。二者均宣称只需提供原始 count 矩阵即可输出 128–768 维的上下文感知嵌入，理论上能显著降低批次效应并提升下游任务精度。

然而，实际落地时却遇到“安装即失败”的连环坑：

硬件门槛：Geneformer 官方建议 40 GB 显存起步（A100-40G），而我本地只有 RTX 3060 12 G；即便开启 gradient_checkpointing + fp16，仍 OOM。

SCGPT 的 PyPI 包依赖 flash-attn>=2.0，该库在 Windows 平台无预编译 whl，源码编译需要 CUDA 11.8+、Visual Studio 2022 与 30 min+ 的 ninja 编译，最终因 cl.exe 与 nvcc 版本冲突而报错退出。

软件生态：Geneformer 的 transformer 版本被锁定在 4.21.1，与 Seurat5 所需的 reticulate 环境里的 4.30+ 冲突；一升级就导致模型权重加载失败（state_dict KeyError）。

SCGPT 的代码仓库最近一次更新已迁移到 Torch 2.1，但配套提供的 scgpt.tokenizer.GeneVocab 类在新版里被重构，旧脚本直接报 AttributeError，而官方文档尚未同步。

 尽管首次尝试以“连包都没装上”告终，但这条思路并不会放弃。后续将继续尝试。

一是轻量化微调：如果全量模型依旧跑不动，就尝试用 LoRA / QLoRA 在 1–2 块 24 GB 卡上做秩为 16 的微调，只训练 adapter 层，预计 6 h 内可收敛。

二是自写“小”模型：若依赖冲突仍无法解决，将基于 transformer encoder 自行实现一个基因排序版 mini-Geneformer（4 层、256 维、4 头），词汇表按 PBMC 高变基因定制，参数仅 5 M，可在 6 G 显存里训练。









