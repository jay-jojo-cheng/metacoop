self_norm_rows <- function(mat){
  num_rows = nrow(mat)
  for (i in 1:num_rows){
    mat[i,] = recenter.rescale(mat[i,])
  }
  return(mat)
}

self_norm_cols <- function(mat){
  num_cols = ncol(mat)
  for (i in 1:num_cols){
    mat[,i] = recenter.rescale(mat[,i])
  }
  return(mat)
}

recenter.rescale <- function(x){
  r.vec = (x - mean(x)) / sd(x)
  return(r.vec)
}
# 
# # alternating as in ichard  Olshen  and  Bala  Rajaratnam.   Successive  normalization  of  rectangular  arrays.Annals of Statistics, 38(3):1638–1664, 2010.
# for(i in 1:100){
#   oldtest = test
#   if ((i%%2) == 0){
#     test = self_norm_cols(test)
#   }
#   else{
#     test = self_norm_rows(test)
#   }
#   print(i)
#   print(sum((test-oldtest)^2))
# }
# 
# # simultaneous as in https://www.jmlr.org/papers/volume16/hastie15a/hastie15a.pdf
# for(i in 1:100){
#   oldtest = test
#   test1 = self_norm_cols(test)
#   test2 = self_norm_rows(test)
#   test = (test1 + test2)/2
#   
#   print(i)
#   print(sum((test-oldtest)^2))
# }
