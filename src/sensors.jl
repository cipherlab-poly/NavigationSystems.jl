"""
    simulateBrownianMotion(PSD, T; dt = T * 1e-3)

returns sample paths of Brownian motion, dB/dt = w, with w white Gaussian
noise with power spectral density psd. T is the time horizon.
dt is the size of the dicrete-time steps used in the simulation.
"""
function simulateBrownianMotion(PSD, T; dt = T*1e-3)
  nsteps = Int(ceil(T/dt))
  res = sqrt(PSD)*randn(nsteps+1)
  res[1] = 0
  return cumsum(res)
end

#b = simulateBrownianMotion(1, 10)
#plot(0:0.01:10,b)
