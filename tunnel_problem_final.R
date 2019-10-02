#Tunnel Problem
tunnel = function() {
  # Set parameters below:
  n = 102 #100 boxes for the tunnel plus 2 for the ends
  # Set final time below and the time step for graphing and saving intermediate values:
  tFinal = 400 #360 to traverse tunnel then dissapation time
  tStep  = 1
  k = 0.01666 * tStep
  va = 0.277 * tStep
  ka = k + va
  
  
  transport = matrix(0, n, n)
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
  #transport[n,n] = -ka
  #transport[1,1] = -k
  
  
  # Solve the problem
  C = matrix(0, n, 1)
  egResult = eigen(transport)
  egValues = egResult$values
  egVectors = egResult$vectors
  egVectorsInv = solve(egResult$vectors)
  
  
  # Set the initial conditions
  #C[2] = 10
  
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
  
  #TESTING WITH EMISSIONS OF PARTICLES CHANGE TO CONCENTRATION LATER
  emission_concentration = 50 / 3600 / (15 * 6 * 100)
  for (j in 1:tSteps) {
    t = tEachStep * j
    timeSeries[[j]] = t
    car_position = t * 1000 / 36
    car_box = floor(car_position / 100) + 2
    #  print(car_box)
    if (j == 1) {
      solutionThisStep = egVectors %*% diag(exp(egValues * t)) %*% egVectorsInv %*% C
      if (car_box <= n) {
        solutionThisStep[car_box] = solutionThisStep[car_box] + emission_concentration *
          tStep
      }
      #     print('yeet')
    }
    
    else{
      solutionThisStep = egVectors %*% diag(exp(egValues) * tStep) %*% egVectorsInv %*% solutionLastStep
      if (car_box <= n) {
        solutionThisStep[car_box] = solutionThisStep[car_box] + emission_concentration *
          tStep
      }
      #    print('doubleyeet')
    }
    solutionThisStep[1] = 0
    solutionThisStep[n] = 0
    #  print(solutionThisStep[car_box])
    #  solutionThisStep[car_box] = solutionThisStep[car_box]+50/3600
    #  print(solutionThisStep[car_box])
    for (k in 1:n) {
      if (j == 1) {
        boxTimeSeries[[k]] = vector("list", tSteps) # create a second dimension
      }
      boxTimeSeries[[k]][[j]] = as.double(solutionThisStep[k]) # remove the imaginary part while we are at it
    }
    timeTimeSeries[[j]] = as.double(solutionThisStep)
    solutionLastStep = timeTimeSeries[[j]]
  }
  
  #PLOTTING
  #boxXaxis = 1:n
  
  # for (j in 1:400) {
  #   plot(
  #     boxXaxis,
  #     timeTimeSeries[[j]],
  #     type = "l",
  #     lwd = 2,
  #     col = "blue",
  #     ylim = c(0, 4*10^(-5)),
  #     main = paste("t =", timeSeries[[j]]),
  #     xlab = "Box#",
  #     ylab = "Concentration"
  #   )
  # }
  #plot(boxXaxis, timeTimeSeries[[tSteps-1]], type="l", lwd=2, col="blue", ylim=c(0, 50), main=paste("t =", timeSeries[[j]]), xlab="Box#", ylab="Concentration")
  return(list(timeTimeSeries, transport, egVectors,egValues))
  #Function returns Green's function, transport matrix, and eigenvalues and eigenvectors of the transport matrix
}