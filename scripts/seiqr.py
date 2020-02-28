from models import seiqr_solver,to_beta
import numpy as np
import sys

"""
Arguments
1. N
2. I0
3. R0
4. t
"""
"""
R0_I = R0_Q
"""

N = int(sys.argv[1])
I0 = int(sys.argv[2])
R0 = float(sys.argv[3])
t = int(sys.argv[4])
dInc = 5.2
dInf = 2.3
epsilon = 0
alpha = 1/dInc
gamma = 1/dInf

"""
y0 = [S0, E0, I0, Q0, R0, CI0, CQ0]
"""
y0 = [N,0,I0,0,0,0,0]

beta_I = to_beta(N=N,durations=dInf,R0=R0)
beta_Q = to_beta(N=N,durations=dInf,R0=R0)

resultDays, resultWk = seiqr_solver(y0=y0,
                                    t=t,
                                    beta_I=beta_I,
                                    beta_Q=beta_Q,
                                    alpha=alpha,
                                    epsilon_E=epsilon,
                                    epsilon_I=epsilon,
                                    gamma=gamma)

print(resultDays.to_dict(orient='Records'))