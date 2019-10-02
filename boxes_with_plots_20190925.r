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
# boxXaxis = 1:n
#matplot(boxTimeSeries[[1]], type="l", lwd=2, ylim=c(0, max(C)), xlab="Time Step", ylab="Concentration")
#legend_text <- c()
#cl<-rainbow(n)
#for (i in 1:n)
#{lines(unlist(boxTimeSeries[[i]]), lwd=2,col=cl[i])
#  legend_text<- c(legend_text,i)
# title(main='Change in Concentration Over Time')
#}
#legend("topright", title="Box Number", legend=legend_text, lwd=c(2,2), col=(rainbow(n)))

# # Plot concentrations for all boxes, x-axis are #
# boxXaxis = 1:n
# for(j in 1:tSteps)
#   {
#   jpeg(paste0(j,".jpeg"))
#   plot(boxXaxis, timeTimeSeries[[j]], type="l", lwd=2, xlab="Box Number" ,ylab = "Concentration", main=paste("t =", timeSeries[[j]]),col="blue", ylim=c(0, max(C)))
#   #paste("myplot_",timeTimeSeries[j], ".jpeg", sep="")
#   #mypath <- file.path("C:/Users/nadiacolombi/Documents/R","R","SAVEHERE",paste("myplot_", timeTimeSeries[[j]], ".jpg", sep = ""))
#   dev.off()}


##   Fit to a Quadratic and solve for D
#   plot(boxXaxis, log(timeTimeSeries[[200]]), type="p", lwd=2,
#   xlab="Box Number", ylab = "log(Concentration)",
#   main=paste("t =", timeSeries[[200]]),col="blue", ylim=c(-20, max(0)))
#   quad <- lm(log(unlist(timeTimeSeries[200]))~poly(boxXaxis,2,raw=TRUE))
#   lines(boxXaxis,predict(quad,data.frame(boxXaxis)),col="red")
#   graphics.off
#
#   poly <- quad[1]$coefficients
#   t <- 200
#   X <- poly[3]
#   #set quadratic equal to f(x,t) to get D
#   D <- (-1/(X*(pi^0.5)*(4*t)^1.5))^(2/3)
#   summary(quad)
#  #Coefs -39.3 + 2.63*X + -0.046*X^2", D=0.006648437
