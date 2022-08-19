# hMKL
we propose a hierarchical multi-kernel learning (hMKL) approach, a novel cancer molecular subtyping method that adopts a two-stage kernel learning strategy to identify cancer subtypes. In stage 1, we obtain a composite kernel borrowing the cancer integration via multi-kernel learning (CIMLR) idea by optimizing the kernel parameters for individual omics data type. In stage 2, we obtain a final fused kernel through a weighted linear combination of individual kernels learned from stage 1 using an unsupervised multiple kernel learning (UMKL) method. Based on the final fusion kernel, k-means clustering is applied to identify cancer subtypes.

# Installation
Before using the hMKL method, several dependent packages need to be downloaded and installed in advance:<br>
install.packages("SIMLR")<br>
install.packages("dplyr")<br>
install.packages("parallel")<br>
install.packages("Matrix")<br>
install.packages("MASS")<br>
install.packages("mixKernel")<br>
install.packages("quadprog")<br>
install.packages("psych")<br>

# Example
This is a basic example which shows you how to solve a common problem:<br>
<br>
data <- exampledata<br>
<br>
```#StageⅠ： Optimize kernel parameters for each omics data type under the CIMLR framework<br>
<br>
mRNA_noc <- CIMLR_noc(data[[1]],cores.ratio = 0)<br>
miRNA_noc <- CIMLR_noc(data[[2]],cores.ratio = 0) <br>
methy_noc <- CIMLR_noc(data[[3]],cores.ratio = 0)<br>
<br>
<br>
Ss_mRNA_noc <- mRNA_noc$S<br>
Ss_miRNA_noc <- miRNA_noc$S<br>
Ss_methy_noc <- methy_noc$S<br>
```
# Stage 2: Obtain the final weighted similarity matrix by UMKL
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

## estimate the optimal number of groups
NUMC = 2:5
res = SIMLR_Estimate_Clusters_W(meta.kernel_sparse_noc$kernel, NUMC = NUMC)
full_K1= NUMC[which.min(res$K1)]
full_K2= NUMC[which.min(res$K2)]

# perform K-means clustering
group_sparse_kmeans_noc = kmeans(meta.kernel_sparse_noc$kernel,full_K1,nstart = 30) 
group_sparse_kmeans <- group_sparse_kmeans_noc$cluster
group_sparse_kmeans
