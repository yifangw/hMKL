#combine.kernels2

combine.kernels2 <- function (..., scale = TRUE, method = c("full-UMKL", "STATIS-UMKL", 
                                        "sparse-UMKL"), knn = 5, rho = 20) 
{
  X <- list(...)
  method <- match.arg(method)
  if (!is.logical(scale)) {
    stop("scale must be either TRUE or FALSE")
  }
  if (!is.list(X)) {
    stop("X must be a list")
  }
  if (length(unique(names(X))) != length(X)) {
    stop("Each block of 'X' must have a unique name.")
  }
  X.classes <- unique(unlist(lapply(X, class)))
  if (length(X.classes) != 1 || X.classes[1] != "kernel") {
    stop("Each block of 'X' must be a kernel.")
  }
  if (length(unique(unlist(lapply(X, function(x) {
    ncol(x$kernel)
  })))) != 1) {
    stop("Unequal number of observations among the kernels of 'X'")
  }
  if (!is.numeric(knn) || knn < 1 || !is.finite(knn)) {
    stop("invalid value for 'knn'.", call. = FALSE)
  }
  if (!is.numeric(rho) || rho < 1 || !is.finite(rho)) {
    stop("invalid value for 'rho'.", call. = FALSE)
  }
  if (scale) {
    X.scaled <- lapply(X, function(x) {
      x.cosinus <- sweep(sweep(x$kernel, 2, sqrt(diag(x$kernel)), 
                               "/"), 1, sqrt(diag(x$kernel)), "/")  #对x$kernel的列执行除以sqrt(diag(x$kernel)的操作，然后再对前面得到的矩阵的行执行除以sqrt(diag(x$kernel)的操作
      t(t(x.cosinus - colSums(x.cosinus)/nrow(x.cosinus)) - 
          rowSums(x.cosinus)/nrow(x.cosinus)) + sum(x.cosinus)/nrow(x.cosinus)^2
    })
  }
  else {
    X.scaled <- list()
    for (i in 1:length(X))
    {
      X.scaled[[i]] <- X[[i]][["kernel"]]
    }
  }
  beta <- 1/length(X.scaled) #初始的beta值为核个数的均值：1/kernels
  if (method == "STATIS-UMKL") {
    similarities <- outer(1:length(X.scaled), 1:length(X.scaled), 
                          FUN = Vectorize(function(i, j) {
                            tr(X.scaled[[i]] %*% X.scaled[[j]])/(norm(X.scaled[[i]], 
                                                                      type = "F") * norm(X.scaled[[j]], type = "F"))
                          }))
    weights <- eigen(similarities, symmetric = TRUE)$vectors[, 
                                                             1]
    weights <- weights/sum(weights)
  }
  else {
############################################################################
#当为sparse-UMKL时的计算
############################################################################    
      all.adjacency <- lapply(X.scaled, function(x) {
      adjacency.x <- apply(x, 1, function(row) {
        adjacency.row <- rank(row)
        adjacency.row[which(adjacency.row < length(adjacency.row) - 
                              knn)] <- 0
        adjacency.row[which(adjacency.row > 0)] <- 1
        adjacency.row
      })
      diag(adjacency.x) <- 0
      adjacency.x + t(adjacency.x)
    })
    graph.weights <- Reduce("+", all.adjacency) #Reduce函数是将每次计算后的结果保留，并与下一个数字进行计算
    C.matrix <- outer(1:length(X.scaled), 1:length(X.scaled),  #outer()求向量的外积
                      FUN = Vectorize(function(i, j) {  #Vectorize()将一个不能进行向量化运算的函数进行转化，使之具备向量化运算功能
                        sum(graph.weights * apply(X.scaled[[i]], 1, function(rowi) {
                          apply(X.scaled[[i]], 1, function(rowj) {
                            sum(abs(rowi - rowj))
                          })
                        }) %*% apply(X.scaled[[j]], 1, function(rowi) { # %*%表示两个矩阵乘积
                          apply(X.scaled[[j]], 1, function(rowj) {
                            sum(abs(rowi - rowj))
                          })
                        }))
                      }))  #得到的矩阵为行数列数均为输入的核矩阵个数的对角阵
    C.matrix <- matrix(unlist(C.matrix), nrow = length(X.scaled))  #行数为核矩阵的个数，列数为样本数，即核矩阵的行列数
    if (method == "sparse-UMKL") {
      dvec <- rep(0, length(X.scaled))
      C.matrix.n <- sweep(sweep(C.matrix, 2, sqrt(diag(C.matrix)), #sweep()：令C.matrix的每一列除以sqrt(diag(C.matrix)，对得到的矩阵的每一行再除以sqrt(diag(C.matrix)
                                "/"), 1, sqrt(diag(C.matrix)), "/")
      Dmat <- 2 * C.matrix.n
      Amat <- matrix(data = 1, nrow = nrow(Dmat), ncol = 1) #生成nrow行、1列，数值均为1的矩阵
      Amat <- cbind(Amat, diag(1, nrow = nrow(Dmat), ncol = nrow(Dmat))) #将行和列均为nrow(Dmat)的对角阵与Amat合并起来
      bvec <- rep(0, length(X.scaled) + 1)
      bvec[1] <- 1
      weights <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)$solution  #求解二次规划问题
    }
    ###################################################################################################
     else {
      k <- 1
      Z <- rep(1/sqrt(length(X.scaled)), length(X.scaled))
      Y <- rep(0, length(X.scaled))
      threshold <- 10^(-2)
      C.matrix.n <- sweep(sweep(C.matrix, 2, sqrt(diag(C.matrix)), 
                                "/"), 1, sqrt(diag(C.matrix)), "/")
      Dmat <- C.matrix.n + diag(rho/2, ncol = length(X.scaled), 
                                nrow = length(X.scaled))
      Dmat <- 2 * Dmat
      Amat <- diag(1, nrow = dim(Dmat)[1], ncol = dim(Dmat)[1])
      bvec <- rep(0, length(X.scaled))
      repeat {
        dvec <- rho * Z - Y
        X.ADMM <- solve.QP(Dmat, dvec, Amat, bvec, meq = 0)$solution
        newZ <- X.ADMM/norm(matrix(X.ADMM), "2")
        if (k != 1 && norm(matrix(Z - newZ)) < threshold) 
          break
        Z <- newZ
        Y <- Y + rho * (X.ADMM - Z)
        k <- k + 1
      }
      weights <- Z/norm(matrix(Z))
    }
  }
  meta.kernel <- lapply(as.list(1:length(X.scaled)), function(x) {
    X.scaled[[x]] * weights[x]
  })
  meta.kernel <- Reduce("+", meta.kernel)
  cl <- match.call()
  cl[[1]] <- as.name("combine.kernels")
  result <- list(call = cl, kernel = meta.kernel, X = X, weights = weights)
  class(result) <- c("kernel", "metaKernel")
  return(invisible(result))
}