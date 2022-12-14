# hMKL
we propose a hierarchical multi-kernel learning (hMKL) approach, a novel cancer molecular subtyping method that adopts a two-stage kernel learning strategy to identify cancer subtypes. In stage 1, we obtain a composite kernel borrowing the cancer integration via multi-kernel learning (CIMLR) idea by optimizing the kernel parameters for individual omics data type. In stage 2, we obtain a final fused kernel through a weighted linear combination of individual kernels learned from stage 1 using an unsupervised multiple kernel learning (UMKL) method. Based on the final fusion kernel, k-means clustering is applied to identify cancer subtypes.

# Installation
Before using the hMKL method, several dependent packages need to be downloaded and installed in advance:<br>
```java
install.packages("SIMLR")
install.packages("dplyr")
install.packages("parallel")
install.packages("Matrix")
install.packages("MASS")
install.packages("mixKernel")
install.packages("quadprog")
install.packages("psych")
```

# Example
This is a basic example which shows you how to solve a common problem:<br>
<br>
```java
data <- exampledata
```

```java

#StageⅠ： Optimize kernel parameters for each omics data type under the CIMLR framework

mRNA_noc <- CIMLR_noc(data[[1]],cores.ratio = 0)

Computing the multiple Kernels.
Performing network diffiusion.
Iteration:  1 
Iteration:  2 
Iteration:  3 
Iteration:  4 
Iteration:  5 
Iteration:  6 
Iteration:  7 
Iteration:  8 
Iteration:  9 
Iteration:  10 
Iteration:  11 
Iteration:  12 
Iteration:  13 
Iteration:  14 
Iteration:  15 
Iteration:  16 
Iteration:  17 
Iteration:  18 
Iteration:  19 
Iteration:  20 
Iteration:  21 
Iteration:  22 
Iteration:  23 
Iteration:  24 
Iteration:  25 
Iteration:  26 
Iteration:  27 
Iteration:  28 
Iteration:  29 
Iteration:  30 

miRNA_noc <- CIMLR_noc(data[[2]],cores.ratio = 0) 

Computing the multiple Kernels.
Performing network diffiusion.
Iteration:  1 
Iteration:  2 
Iteration:  3 
Iteration:  4 
Iteration:  5 
Iteration:  6 
Iteration:  7 
Iteration:  8 
Iteration:  9 
Iteration:  10 
Iteration:  11 
Iteration:  12 
Iteration:  13 
Iteration:  14 
Iteration:  15 
Iteration:  16 
Iteration:  17 
Iteration:  18 
Iteration:  19 
Iteration:  20 
Iteration:  21 
Iteration:  22 
Iteration:  23 
Iteration:  24 
Iteration:  25 
Iteration:  26 
Iteration:  27 
Iteration:  28 
Iteration:  29 
Iteration:  30 

methy_noc <- CIMLR_noc(data[[3]],cores.ratio = 0)

Computing the multiple Kernels.
Performing network diffiusion.
Iteration:  1 
Iteration:  2 
Iteration:  3 
Iteration:  4 
Iteration:  5 
Iteration:  6 
Iteration:  7 
Iteration:  8 
Iteration:  9 
Iteration:  10 
Iteration:  11 
Iteration:  12 
Iteration:  13 
Iteration:  14 
Iteration:  15 
Iteration:  16 
Iteration:  17 
Iteration:  18 
Iteration:  19 
Iteration:  20 
Iteration:  21 
Iteration:  22 
Iteration:  23 
Iteration:  24 
Iteration:  25 
Iteration:  26 
Iteration:  27 
Iteration:  28 
Iteration:  29 
Iteration:  30 

Ss_mRNA_noc <- mRNA_noc$S
Ss_miRNA_noc <- miRNA_noc$S
Ss_methy_noc <- methy_noc$S
```
```java
#Stage 2: Obtain the final weighted similarity matrix by UMKL
set_input <- function(kernel_mat)
{
  output <- list(kernel = kernel_mat)
  class(output) <- "kernel"
  return(output)
}



mRNA_kernel_noc <- set_input(Ss_mRNA_noc)
miRNA_kernel_noc <- set_input(Ss_miRNA_noc)
methy_kernel_noc <- set_input(Ss_methy_noc)

meta.kernel_sparse_noc <- combine.kernels2(x = mRNA_kernel_noc,y = miRNA_kernel_noc,z = methy_kernel_noc,
                                           method = "sparse-UMKL", scale=F)
```
```java
# estimate the optimal number of groups
NUMC = 2:5
res = SIMLR_Estimate_Clusters_W(meta.kernel_sparse_noc$kernel, NUMC = NUMC)
full_K1= NUMC[which.min(res$K1)]
full_K2= NUMC[which.min(res$K2)]
```
```java
# perform K-means clustering
group_sparse_kmeans_noc = kmeans(meta.kernel_sparse_noc$kernel,full_K1,nstart = 30) 
group_sparse_kmeans <- group_sparse_kmeans_noc$cluster
group_sparse_kmeans
```
