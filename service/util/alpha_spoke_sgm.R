library(igraph)

# Simulates an nxn permutation matrix uniformly
# INPUT: n : dimension of the desired matrix
# OUTPUT: nxn permutation matrix
gen_rand_perm_mat <- function(n){
    perm <- sample(n)
    mat <- Matrix::Diagonal(n)
    mat <- mat[perm,]
    return(mat)
}

# Simulate an alpha-spoke doubly stochastic matrix using the following procedure:
# (1) uniformly generate a permutation matrix P
# (2) Let the barycenter doubly stochastic matrix be B, return the matrix
#     (1-alpha)*B + alpha*P
# INPUT: n : dimension of desired matrix
#        alpha : distance along line segment between barycenter and permutation matrix
# OUTPUT: nxn doubly stochastic matrix
gen_rand_ds_mat <- function(n, alpha){
        P <- gen_rand_perm_mat(n)
        B <- Matrix::Matrix(1/n, nrow=n, ncol=n)
        return((1 - alpha)*B + alpha*P)
}

# One iteration of the alpha-spoke random initialization
# INPUT: A,B : adjacency matrices of two graphs with equal number of vertices
#        m : number of seeds
#        alpha : relative distance of initialization from the barycenter
#        num_starts : number of runs of SGM to be used in the average
#        fw_iter : max. number of Frank-Wolfe iterations used in FAQ
# OUTPUT: A doubly stochastic matrix where the i,j-th entry represents the
#         probability that vertex i in graph A is matched to vertex j in
#         graph B.

one_rand_start <- function(A, B, m, alpha, fw_iter){
    n <- dim(A)[1]
    start <- as.matrix(gen_rand_ds_mat(n-m, alpha))
    P <- match_vertices(A,B, m, start=start, fw_iter)
    return(P$P)
}


# Runs SGM with a random alpha-spoke matrix as initialization and returns the average
# of the permutation matrices returned by SGM.
# INPUT: A,B : adjacency matrices of two graphs with equal number of vertices
#        m : number of seeds
#        alpha : relative distance of initialization from the barycenter
#        num_starts : number of runs of SGM to be used in the average
#        fw_iter : max. number of Frank-Wolfe iterations used in FAQ
# OUTPUT: A doubly stochastic matrix where the i,j-th entry represents the
#         probability that vertex i in graph A is matched to vertex j in
#         graph B.
mv_rand_start <- function(A,B, m, alpha, num_starts, fw_iter){
        res.mat <- foreach(i=1:num_starts, .combine='+') %do% one_rand_start(A, B, m, alpha, fw_iter)
        return(res.mat/num_starts)
}
####################################################################################
#                               FOR TESTING                                        #
####################################################################################
#set.seed(144169)
#g1 <- erdos.renyi.game(20, .2)
#randperm <- c(1:3, 3 + sample(17))
#g2 <- sample_correlated_gnp(g1, corr=0.8, p=g1$p, perm=randperm)
#A <- as.matrix(get.adjacency(g1))
#B <- as.matrix(get.adjacency(g2))

#P <- mv_rand_start(A, B, 3, 0.1, 10, 10)


compare_rand_start <- function(n, p, rho, m, alpha=0.1, num_starts=100, fw_iter=20, myseed=-1){
        if (myseed != -1){
                set.seed(myseed)
        }
        g1 <- erdos.renyi.game(n, p)
        randperm <- c(1:m, m + sample(n-m))
        g2 <- sample_correlated_gnp(g1, corr=rho, p=g1$p, perm=randperm)
        A <- as.matrix(get.adjacency(g1))
        B <- as.matrix(get.adjacency(g2))

        barycenter <- matrix(1/(n-m), nrow=n-m, ncol=n-m)
        start.time <- proc.time()
        res.sgm <- match_vertices(A, B, m, start=barycenter, fw_iter)
        end.sgm <- proc.time()
        res.rand.start <- mv_rand_start(A, B, m, alpha, num_starts, fw_iter)
        end.rand.starts <- proc.time()
        #print("Time elapsed for SGM:\n")
        #print(end.sgm - start.time)
        #print("Time elapsed for random start SGM:\n")
        #print(end.rand.starts - end.sgm)

        #print("Number of errors for SGM: ")
        #print(n-m - sum(res.sgm$corr[,2] == randperm[-1:-m]))
        #print("Number of errors for random start SGM: ")
        guess.perm <- apply(res.rand.start, 1, which.max)
        #print(n-m - sum(guess.perm[-1:-m] == randperm[-1:-m]))
        sgm.num.correct <- sum(res.sgm$corr[,2] == randperm[-1:-m])
        rand.start.num.correct <- sum(guess.perm[-1:-m] == randperm[-1:-m])
        return(list(SGM=sgm.num.correct/(n-m), RSSGM=rand.start.num.correct/(n-m)))
}
count_correct_sgm <- function(res.sgm, randperm, m){
        return(sum(res.sgm$corr[,2] == randperm[-1:-m]))
}
# k = number of potential vertex matches to look at
count_correct_rssgm <- function(res.rssgm, randperm, m, k){
        #guess.perm <- apply(res.rssgm, 1, which.max)
        num.vertices <- dim(res.rssgm)[1]
        potential.matches <- apply(res.rssgm[(m+1):num.vertices,(m+1):num.vertices], 1, order)[(num.vertices-m):(num.vertices-m-k+1),]
        potential.matches <- as.vector(t(potential.matches))
        return(sum((potential.matches+m) == randperm[-1:-m]))
}

sample_corr_erdos_renyi <- function(n, p, rho, m, myseed=-1){
        if (myseed != -1){
                set.seed(myseed)
        }
        g1 <- erdos.renyi.game(n, p)
        randperm <- c(1:m, m + sample(n-m))
        g2 <- sample_correlated_gnp(g1, corr=rho, p=g1$p, perm=randperm)
        A <- as.matrix(get.adjacency(g1))
        B <- as.matrix(get.adjacency(g2))
        return(list(A=A, B=B, permutation=randperm))
}

gen_barycenter <- function(k){
        return(matrix(1/k, nrow=k, ncol=k))
}


run_simulation <- function(n, p, rho, m, alpha=0.1, num_starts=20, fw_iter=20, num_sim=10){
        sum.sgm.correct <- 0
        sum.rssgm.correct <- 0
        for(k in 1:num_sim){
                res <- compare_rand_start(n, p, rho, m, alpha, num_starts, fw_iter)
                sum.sgm.correct = sum.sgm.correct + res$SGM
                sum.rssgm.correct = sum.rssgm.correct + res$RSSGM
        }
        return(list(SGM=sum.sgm.correct/num_sim, RSSGM=sum.rssgm.correct/num_sim))
}
big_sim <- function(n=100, p=0.1, rho=c(0.2, 0.5), num.seeds=c(10,20), alpha=c(0.1, 0.3, 0.5), num.starts=c(10,20), num.sim=10, fw_iter=20){
        function.call <- match.call

        R <- length(rho)
        S <- length(num.seeds)
        AL <- length(alpha)
        N <- length(num.starts)
        sgm.score <- array(0, dim=c(R, S, AL, N), dimnames=c("rho", "m", "alpha", "num.starts"))

        rssgm.score.1 <- array(0, dim=c(R, S, AL, N), dimnames=c("rho", "m", "alpha", "num.starts"))
        rssgm.score.2 <- array(0, dim=c(R, S, AL, N), dimnames=c("rho", "m", "alpha", "num.starts"))
        for(sim in 1:num.sim){
                print(paste("Simulation #", sim))
        for(i in 1:R){
        for(j in 1:S){
                corr.ER <- sample_corr_erdos_renyi(n, p, rho[i], num.seeds[j])
                A <- corr.ER$A
                B <- corr.ER$B
                res.sgm <- match_vertices(A, B, num.seeds[j],
                                          start=gen_barycenter(n-num.seeds[j]),
                                          iteration=fw_iter)
                sgm.correct <- count_correct_sgm(res.sgm,
                                              corr.ER$permutation,
                                              num.seeds[j])
                for(k in 1:AL){
                for(l in 1:N){
                        sgm.score[i,j,k,l] <- sgm.score[i,j,k,l]+ sgm.correct/(n-num.seeds[j])
                                                


                        time.start <- proc.time()
                        res.rssgm <- mv_rand_start(A, B, num.seeds[j],
                                                  alpha[k],
                                                 num.starts[l],
                                                fw_iter) 
                        rssgm.correct <- count_correct_rssgm(res.rssgm,
                                                          corr.ER$permutation,
                                                          num.seeds[j],
                                                          1)
                        rssgm.score.1[i,j,k,l] <- rssgm.score.1[i,j,k,l] + rssgm.correct/(n-num.seeds[j])
                        rssgm.correct <- count_correct_rssgm(res.rssgm,
                                                          corr.ER$permutation,
                                                          num.seeds[j],
                                                          2)

                        rssgm.score.2[i,j,k,l] <- rssgm.score.1[i,j,k,l] + rssgm.correct/(n-num.seeds[j])                      
                        time.end <- proc.time()
                        out.str <- paste("rho=",rho[i],
                                       "number of seeds=", num.seeds[j],
                                       "alpha=", alpha[k],
                                       "number of starts=", num.starts[l])
                        print(out.str)
                        print(time.end-time.start)
                }
                }
        }
        }
        }
        sgm.score <- sgm.score/num.sim
        rssgm.score.1 <- rssgm.score.1/num.sim
        rssgm.score.2 <- rssgm.score.2/num.sim
        return(list(call=function.call, 
                    SGM=sgm.score, 
                    RSSGM.1=rssgm.score.1, 
                    RSSGM.2=rssgm.score.2))
}

plot_sgm <- function(x, y1, y2, col1, col2, main.title, label){
        min.val <- min(c(y1,y2))
        max.val <- max(c(y1,y2))
        plot(x,y1, ylim=c(min.val, max.val), col=col1, type="b", pch=16,  
             main=main.title,
             xlab=label,
             ylab="Percent Correctly Matched Vertices")
        lines(x,y2, col=col2, type="b", pch=16)
        legend("topright", col=c(col1, col2), legend=c("sgm", "rssgm"),
               pch=16)
}
plot_rssgm <- function(res, rho, num.seeds, alpha, num.init){
        R <- length(rho)
        N <- length(num.seeds)
        A <- length(alpha)
        I <- length(num.init)

        par(mfrow=c(R,N), mar=c(5,5,1,1), oma=c(1,1,3,1))

        for(r in 1:R){
        for(n in 1:N){
                title <- paste("rho=", rho[r], ", #seeds=", num.seeds[n]) 
                m <- min(c(res$SGM[r,n,,1], res$RSSGM.1[r,n,,1], res$RSSGM.2[r,n,,1]))
                M <- max(c(res$SGM[r,n,,1], res$RSSGM.1[r,n,,1], res$RSSGM.2[r,n,,1]))
                plot(alpha, res$SGM[r,n,,1], ylim=c(m,M), col="red", type="b", pch=16,
                     main=title,
                     xlab="alpha",
                     ylab="Percent Correctly Matched Vertices")
                lines(alpha, res$RSSGM.1[r,n,,1], col="blue", type="b", pch=16)
                lines(alpha, res$RSSGM.2[r,n,,1], col="green", type="b", pch=16)
                legend("bottomleft", col=c("red","blue", "green"),
                       legend=c("sgm", "rssgm 1", "rssgm 2"),
                       pch=16)
                mtext("n=100, p=0.1, #restarts=20, #simulations=20, fw_iter=10",
                      side=3, line=1, outer=TRUE)
        }
        }
}

