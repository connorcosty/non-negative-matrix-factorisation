rm(list = ls())
# Load the data, make sure that you are in the right folder
load(file = "CW2_25.RData")
uob_email <- "AB12345@bristol.ac.uk" 
#university email obscured
my_seed <- as.numeric(gsub("\\D", "", uob_email))
print(my_seed)
set.seed(my_seed)
A_test <- matrix(sample(1:9), 3, 3)
V_test <- matrix(sample(1:12), 4, 3)
W_test <- matrix(sample(1:8), 4, 2)
H_test <- matrix(sample(1:6), 2, 3)
my_seed_2 <- sample(c(1, 3, 4, 5, 6))[1]
set.seed(my_seed_2)
print(my_seed_2)
V <- V[sample(nrow(V)), ]
class(V)
dim(V)
range(V)
print(V)
#The above code loads the data from the file into V, a 400x4096 data matrix.
imageplot <- function(X, title = ''){
  Xrev <- t(apply(X, 2, rev))
  image(Xrev, col = gray(seq(0,1,length=256)), main = title)
}
#Image plotting function requires a 1 dimensional input for X

x <- 1:9 / 9
X <- matrix(x,nrow = 3, ncol = 3)
imageplot(X, title = "Test Image")

plot_image <- function(x, title){
  imageplot(matrix(V[x, ],nrow=64,ncol=64), title)
}
plot_image(1, "row 1")
#The above function plots a 64x64 image from a row of the V data matrix

mat_init <- function(n, d, r, u1, u2){
  W <- matrix(u1[1:(n*r)], nrow = n, ncol = r)
  H <- matrix(u2[1:(d*r)], nrow = r, ncol = d)
  return(list(W,H))
}
#The above function initialises W and H to equal input vectors u1 and u2.

set.seed(123)
u1 <- runif(1e6)
u2 <- runif(1e6)
mat_init(n=4, d=3, r=2, u1 = u1, u2 = u2)
#The above code will set the seed for random number generation to 123,
#then will set u1 and u2 to a uniform distribution given by the seed
#then will initialise W and H

frob <- function(A){
  return(sqrt(sum(A*A)))
}
frob(A = A_test)
#The above function returns the Frobenius Norm of an input matrix

err_frob <- function(V, W, H){
  return(frob(V- (W%*%H)))
}
#The above function calculates the Frobenius norm of the error matrix
#where the error matrix is the matrix of element wise differences 
#between V and the W*H approximation
err_frob(V = V_test, W = W_test, H = H_test)

up_H <- function(V, W, H, eps){
  H <- H * ( (t(W)%*%V) / (((t(W) %*%W) %*%H) + eps ))
}
up_W <- function(V, W, H, eps){
  W <- W * ( (V%*% t(H)) / (((W%*%H) %*%t(H)) + eps) )
}
#The above two functions will update the approximation of W and H
up_H(V = V_test, W = W_test, H = H_test, eps = 1e-9)
up_W(V = V_test, W = W_test, H = H_test, eps = 1e-9)

lowrank_fit <- function(V, r, u1, u2, maxit, delta, eps){
  initial <- mat_init(400, 4096, r, u1, u2)
  W <- initial[[1]]
  H <- initial[[2]]
  norms <- list()
  alpha <- Inf
  #The above code initialises W, H, norms and alpha
  while((alpha - err_frob(V,W,H))/err_frob(V,W,H) > delta){
    #The above calculates if the change of the Frobenius norm of the
    #Error matrix is sufficiently close to delta, exiting the function
    #when it is sufficiently close
    alpha <- err_frob(V,W,H)
    #Alpha is the previous Frobenius norm of the Error matrix
    H <- up_H(V,W,H,eps)
    W <- up_W(V,W,H,eps)
    #Updates H and W
    norms <- c(norms, err_frob(V,W,H))
    #Appends the next Frobenius norm of the Error matrix to the 
    #list 'norms'
    if(length(norms) %% 10 == 0){
      cat(".")
      #Prints . every 10th operation
      if(length(norms) %% 100 == 0){
        cat("\n")
        #Goes to a new line every 100th operation
      }
    }
    if(maxit == length(norms)){
      print("Maxit reached!")
      break
      #breaks the loop if the maximum iterations is reached
      #prints as such
    }
  }
  if(maxit != length(norms)){
    print("Convergence")
    #If the loops ends without maxit = length(norms) then
    #the function ended prematurely, so converges
  }
  return(list(list(W,H), norms))
  #returns W and H in a single object, then the list of Frobenius
  #norms of error matrices as a single object
}

tic <- proc.time()
fit16 <- lowrank_fit(V = V, r = 4, u1 = u1, u2 = u2, 
                     maxit = 500, delta = 1e-3, eps = 1e-9)
proc.time() - tic
#Lowrank_fit for r=4 (64/4=16 hence the name fit16, however the
#naming system R2, R8, R16 is easier for me to understand hence its use)

Iterations <- 1:length(fit16[[2]])
ReconstructionError <- fit16[[2]]
plot(Iterations, ReconstructionError)

#near 129 the function breaks prematurely, so converges at about that point,
#implying that the reconstruction error won't reduce significantly with more 
#iterations

R2 <- lowrank_fit(V=V, r=2, u1=u1, u2=u2, maxit=500, delta = 1e-3, eps=1e-9)
R8 <- lowrank_fit(V=V, r=8, u1=u1, u2=u2, maxit=500, delta = 1e-3, eps=1e-9)
R16 <- lowrank_fit(V=V, r=16, u1=u1, u2=u2, maxit=500, delta = 1e-3, eps=1e-9)
#the above runs the lowrank_fit function for r=2, 8, and 16, assigning
#the outputs to respectively named variables

plot(c(2, 4, 8, 16), c(R2[[2]][length(R2[[2]])], fit16[[2]][length(fit16[[2]])], 
                       R8[[2]][length(R8[[2]])], R16[[2]][length(R16[[2]])]))
#The code above plots the final frobenius norm of each lowrank_fit output
#against its respective R value.

reconstructed2 <- R2[[1]][[1]] %*% R2[[1]][[2]]
reconstructed4 <- fit16[[1]][[1]] %*% fit16[[1]][[2]]
reconstructed8 <- R8[[1]][[1]] %*% R8[[1]][[2]]
reconstructed16 <- R16[[1]][[1]] %*% R16[[1]][[2]]
#The above assigns the reconstructed image matrix to an appropriately
#named variable

plot.new()
#empties the plot viewer
par(mfrow = c(3,5))
#sets the dimensions of the plot viewer
for(i in 1:3){
  plot_image(i, paste("Original Image", i))
  imageplot(matrix(reconstructed16[i, ],nrow=64,ncol=64), "Reconstructed, R=16")
  imageplot(matrix(reconstructed8[i, ],nrow=64,ncol=64), "Reconstructed, R=8")
  imageplot(matrix(reconstructed4[i, ],nrow=64,ncol=64), "Reconstructed, R=4")
  imageplot(matrix(reconstructed2[i, ],nrow=64,ncol=64), "Reconstructed, R=2")
}
#The above iterates through plotting the original image stored in V (
#hence the plot_image function), then decreasing quality image matrices
#as the value of R plotted decreases from 16 to 2, before the next image
#is then plotted

#The results are interesting; while the increased compression leads to
#lower quality images on the output, despite being unidentifiable,
#each image matrix is still a human face, even when only a 32nd of
#the data is stored. Overall, the trend clearly emerges where the 
#lower the r value, the lower quality is retained in the restored image.
