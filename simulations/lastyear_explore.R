# exploration from last year

# look at treatment,  Y and Y given A=1, Y given  A=0

# getting names in matrix to call in and store y and z values 

addon = rbind(c(0,0,0),c(0,0,1),c(0,1,0),c(0,1,1),c(1,0,0),c(1,0,1),c(1,1,0),c(1,1,1))
addon
dataset_names = do.call(cbind,lapply(1:4,FUN = function(i) {
  vapply(1:8, FUN=function(nm) {
    y=paste0("inbound/pre_data/",paste0(i,"."),paste(addon[nm,],collapse=""),".y.csv")
    z=paste0("inbound/pre_data/",paste0(i,"."),paste(addon[nm,],collapse=""),".z.csv")
    return(c(z,y))
  },FUN.VALUE = c("adsf","sdf"))
}))

# list of dataframes per run
AYlist = lapply(1:32,FUN=function(x) {
  A=read.csv(dataset_names[1,x])
  Y=read.csv(dataset_names[2,x])
  return(data.frame(A=A,Y=Y))
})

# get basic info
basics = data.frame(t(sapply(AYlist, FUN = function(x) {
  Abar = mean(x$z)
  Ybar = mean(x$y)
  Ysd = sd(x$y)
  Ygiven1 = sum(with(x,z*y))/sum(x$z)
  Ygiven0 = sum(with(x,(1-z)*y))/sum(1-x$z)
  a = min(x$y)
  b = max(x$y)
  ATE = Ygiven1-Ygiven0
  return(c(Abar,Ybar,Ysd,Ygiven1,Ygiven0,ATE,a,b))
})))

colnames(basics) = c("Abar","Ybar","Ysd","Ygiven1","Ygiven0","ATE","minY","maxY")

basics

# look at summary of covariates for both y and z
Xy = read.csv("inbound/pre_data/X_subset_y.csv")
Xz = read.csv("inbound/pre_data/X_subset_z.csv")
head(Xy)
ncol(Xy)
summary(Xy)

head(Xz)
ncol(Xz)
summary(Xz)


