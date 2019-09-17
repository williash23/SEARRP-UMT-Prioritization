
# Model to calculate metapopulation capacity for a given landscape matrix



rm(list=ls(all=TRUE))


#setwd("C:/Users/JB/Desktop/Metapop capacity")
setwd("C:/Users/saraw/Documents/Prioritization/Metapop_Capacity")



#------------------- LOAD LANDSCAPE MATRICES -------------------------------
list.files()

# Connectivity matrix
#load("prob_mat_conn_prior.Rdata")
#print(load("prob_mat_conn_prior.Rdata"))
#conn <- prob_mat_conn_prior

load("prob_mat_conn_prior_butterflies.Rdata")
print(load("prob_mat_conn_prior_butterflies.Rdata"))
conn <- prob_mat_conn_prior


# No connectivity matrix
#load("prob_mat_no_conn_prior.Rdata")
#print(load("prob_mat_no_conn_prior.Rdata"))
#noconn <- prob_mat_no_conn_prior

load("prob_mat_no_conn_prior_butterflies.Rdata")
print(load("prob_mat_no_conn_prior_butterflies.Rdata"))
noconn <- prob_mat_no_conn_prior




#------------------- CHECK IF THE TWO MATRICES ARE DIFFERENT ---------------
matequal <- function(x, y)
	is.matrix(x) && is.matrix(y) && dim(x) == dim(y) && all(x == y)
matequal(conn, noconn)  # should be FALSE




#------------------- METAPOPULATION CAPACITY -------------------------------
# Connectivity
mat1 <- conn
patch.areas1 <- diag(mat1)
diag(mat1) <- 1
mattmp <- (mat1) * (patch.areas1) %*% t(patch.areas1)^0.5
MCconn <- Re(eigen(mat1)$values[1]) # Schnell method
#MCconn <- eigen(conn)$values[1] # Hanski method


# No connectivity
mat2 <- noconn
patch.areas2 <- diag(mat2)
diag(mat2) <- 1
mattmp <- (mat2) * (patch.areas2) %*% t(patch.areas2)^0.5
MCnoconn <- Re(eigen(mat2)$values[1]) # Schnell method
#MCnoconn <- eigen(noconn)$values[1] # Hanski method


# Proportional difference
MCconn / MCnoconn







#------------------- CUMULATIVE DISPERSAL -------------------------------
# This is a simple metric of how much dispersal can occur across the whole landscape
# It's just the sum of all the dispersal elements in the landscape matrix


# Connectivity
diag(mat1) <- 0
cs1 <- colSums(mat1)
SumConn <- sum(cs1)


# No connectivity
diag(mat2) <- 0
cs2 <- colSums(mat2)
SumNoConn <- sum(cs2)


# Proportional difference
SumConn / SumNoConn























