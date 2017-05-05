# exploration from last year

# look at treatment,  Y and Y given A=1, Y given  A=0

# getting names in matrix to call in and store y and z values 
UU = rmultinom(10000,3,c(.3,.4,.3))
addon = unique(apply(UU,1,FUN =function(x) {
  as.numeric(x==0)}))
addon = rbind(addon, c(0,0,0))
dataset_names = do.call(cbind,lapply(1:4,FUN = function(i) {
  vapply(1:8, FUN=function(nm) {
    y=paste0("inbound/pre_data/",paste0(c(i),"."),paste(addon[nm,],collapse=""),".y.csv")
    z=paste0("inbound/pre_data/",paste0(c(i),"."),paste(addon[nm,],collapse=""),".z.csv")
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
  ATE = Ygiven1-Ygiven0
  return(c(Abar,Ybar,Ysd,Ygiven1,Ygiven0,ATE))
})))

colnames(basics) = c("Abar","Ybar","Ysd","Ygiven1","Ygiven0","ATE")

basics