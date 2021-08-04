### Create dummy data with 15 instances for SM1, SM2, and DM
### Generate the data from a random uniform distribution bounded by 0 and 1
data <- data.frame(SM1 = runif(15, 0, 1),
                   SM2 = runif(15, 0, 1),
                   DM = runif(15, 0, 1))
### Create a vector of operations to apply to the data frame
OP.VECT <- c()
i <- 1
for (OP in c("+", "-")){
  OP.VECT[i] <- paste("SM1", OP, "SM2")
  i <- i + 1
}

