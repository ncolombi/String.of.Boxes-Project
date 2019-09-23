k = 0.01
va = 0.004
ka = k+va

transport = matrix(c(-ka-k,k,0,ka,ka,-ka-k,k,0,0,ka,-ka-k,k,k,0,ka,-ka-k), nrow=4,ncol=4)
#now let's do n boxes
k = 0.01 #arbitrary
va = 0.014 #arbitrary
ka = k+va
n = 5#insert number of boxes here
transport = matrix(0,n,n)
for (i in 1:n) {
  transport[i,i] = -ka-k
  if(i+1<n){
  transport[i,i+1] = k}
  else{transport[i,1] = k}
  if(i-1>0){
    transport[i,i-1] = ka}
  else{transport[i,n] = ka}
}
C = matrix(0,1,n)
transport
eigen(transport)

