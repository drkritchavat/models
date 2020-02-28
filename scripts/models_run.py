from models import seiqr_solver,to_beta
import numpy as np
import sys
import json

model = sys.argv[1]

filename = f'../params/params.json'
with open(filename,'r') as f:
    params = json.loads(f.read())
for param in params:
    if param['model'] == model:
        selected_model = param
print(selected_model)
N = selected_model['N']
I0 = selected_model['I0']
R0_I = selected_model['R0_I']
R0_Q = selected_model['R0_Q']
t = selected_model['t']
dInc = 5.2
dInf = 2.3
epsilon = 0
alpha = 1/dInc
gamma = 1/dInf

"""
y0 = [S0, E0, I0, Q0, R0, CI0, CQ0]
"""
y0 = [N,0,I0,0,0,0,0]
print(y0)
beta_I = to_beta(N=N,durations=dInf,R0=R0_I)
beta_Q = to_beta(N=N,durations=dInf,R0=R0_Q)

resultDays, resultWk = seiqr_solver(y0=y0,
                                    t=t,
                                    beta_I=beta_I,
                                    beta_Q=beta_Q,
                                    alpha=alpha,
                                    epsilon_E=epsilon,
                                    epsilon_I=epsilon,
                                    gamma=gamma)

print(resultDays[['CI']].to_dict(orient='Records'))
