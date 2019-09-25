# Set parameters below:
k = 0.01
va = 0.014
ka = k+va
n = 25
# Set final time below and the time step for graphing and saving intermediate values:
tFinal = 200
tStep  = 20


transport = matrix(0,n,n)
for (j in 1:n) {
  transport[j, j] = -ka - k
  if (j + 1 <= n) {
    transport[j, j + 1] = k
  }
  else{
    transport[j, 1] = k
  }
  if (j - 1 > 0) {
    transport[j, j - 1] = ka
  }
  else{
    transport[j, n] = ka
  }
}

# Solve the problem
C = matrix(0,n,1)
egResult = eigen(transport)
egValues = egResult$values
egVectors = egResult$vectors
egVectorsInv = solve(egResult$vectors)

# Set the initial conditions
C[2] = 10

# Generate the time series for graphing. 
tSteps = as.integer(tFinal / tStep)
tEachStep = tFinal / tSteps

# Each box's time series will be stored in boxTimeSeries[boxNumber][timeStep]
# An alternative data structure is         timeTimeSeries[timeStep][boxNumber]
# This can be used depending on which kind of plot you want...
# And the steps themselves are in timeSeries[timeStep]
timeSeries = vector("list", tSteps)
timeTimeSeries = vector("list", tSteps)
boxTimeSeries = vector("list", n)


for(j in 1:tSteps) {
  t = tEachStep * j
  timeSeries[[j]] = t
  solutionThisStep = egVectors %*% diag(exp(egValues * t)) %*% egVectorsInv %*% C
  
  for(k in 1:n) {
    if(j == 0) {
      boxTimeSeries[[k]] = vector("list", tSteps) # create a second dimension
    }
    boxTimeSeries[[k]][[j]] = as.double(solutionThisStep[k]) # remove the imaginary part while we are at it
  }
  timeTimeSeries[[j]] = as.double(solutionThisStep)
}

# Plot for each box...
# plot(timeTimeSeries, boxTimeSeries[[1]], type="l", lwd=2, col="blue", ylim=c(0, max(C)), xaxs="i", yaxs="i")
# lines(timeTimeSeries, boxTimeSeries[[2]], lwd=2, col="red")
# lines(timeTimeSeries, boxTimeSeries[[3]], lwd=2, col="green")
# legend("topleft", legend=c("box 1","box 2","box 3"), lwd=c(2,2), col=c("blue","red","green"))

# Plot concentrations for all boxes, x-axis are #
boxXaxis = 1:n

for(j in 1:tSteps) {
  plot(boxXaxis, timeTimeSeries[[j]], type="l", lwd=2, col="blue", ylim=c(0, max(C)), main=paste("t =", timeSeries[[j]]), xlab="Box#", ylab="Concentration")
}
