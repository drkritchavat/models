"""
THAILAND
Simple SEIR model
Density dependence
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def to_beta(R0,durations,N):
    """
    R0: Basic Reproductive Number
    Return beta: per capita rate at which individuals make effective contact 
    (Probability that individuals make effective contact that result in transmission per time unit)
    """
    return (R0/N)/durations

def seiqr_equations(y0,
                    t,
                    beta_I,
                    beta_Q,
                    alpha,
                    epsilon_E,
                    epsilon_I,
                    gamma):
    """
    Modified version of SEIQR eqations (include Q compartment)
    """
    S, E, I, Q, R, CI, CQ = y0
    dS = -beta_I*S*I - beta_Q*S*Q
    dE = beta_I*S*I + beta_Q*S*Q - alpha*E
    dI = (1-epsilon_E)*alpha*E - epsilon_I*I - gamma*I
    dQ = epsilon_E*alpha*E + epsilon_I*I - gamma*Q
    dR = gamma*(I+Q)
    dCI = alpha*E
    dCQ = epsilon_E*E + epsilon_I*I
    return dS, dE, dI, dQ, dR, dCI, dCQ


def seiqr_solver(y0,
                 t,
                 beta_I,
                 beta_Q,
                 alpha,
                 epsilon_E,
                 epsilon_I,
                 gamma):
    """
    SEIQR Solver
    """
    args = (beta_I,
            beta_Q,
            alpha,
            epsilon_E,
            epsilon_I,
            gamma)

    t = np.arange(0,t,1)
    df = pd.DataFrame(
        odeint(seiqr_equations,y0,t,args))

    df = df.reset_index()
    df.columns = ['DAYS','S','E','I','Q','R','CI','CQ']
    dCICQ = df.loc[:,['CI','CQ']].diff()
    dCICQ.columns = ['dCI','dCQ']
    df = pd.concat([df,dCICQ],axis=1)

    epicurve = df.loc[:,['dCI']]
    epicurve['week'] = np.floor(epicurve.index/7)
    epicurve = epicurve.groupby('week')['dCI'].sum()
    
    return df,epicurve


def thai_seir(Iu0,Ic0,Cum0,Cum_local0,ndays,beta,beta2,gamma,N,alpha,impE,impI,pImpc,pIc,pIu_Ic,dr,cr):
    
    """
    C: The number of cumulative infectious cases
    impE: The number of imported latent infection per day
    impI: The number of imported infectious cases per day
    """
    
    init = {'S':N,'E':0,'Iu':Iu0,'Ic':Ic0,'R1':0,'R2':0,'D':0,'C':0,'Cum_det':0,'Cum':Cum0,'Cum_local':Cum_local0}
    t = np.linspace(0, ndays, ndays)
    
    """
    The SEIR model differential equations.
    """
     
    def deriv(y, t, N, beta, beta2, gamma, alpha, impE,impI,pImpc,pIc,pIu_Ic,dr,cr):
       

        S,E,Iu,Ic,R1,R2,C,D,Cum_det,Cum,Cum_local = y
        dSdt = - (beta * S * Iu) - (beta2 * S * Ic)
        dEdt = (beta * S * Iu) + (beta2 * S * Ic) - (alpha * E) + impE
        dIudt = alpha * E * (1-pIc)  + impI * (1-pImpc) - pIu_Ic * Iu - gamma * Iu
        dIcdt = (alpha * E * pIc)  + (impI * pImpc) + (pIu_Ic * Iu) - (gamma * Ic)
        dR1dt = gamma * Iu
        dR2dt = gamma * Ic - dr * R2 - cr * R2
        dDdt = dr * R2
        dCdt = cr * R2
        dCum_det_dt = (alpha * E * pIc)  + (impI * pImpc) + (pIu_Ic * Iu) # Cumulative reported cases (which has been detected and isolated)
        dCumdt = (alpha * E) + impI # Cumulative total infectious cases
        dCum_localdt = beta * S * Iu # Cumulative total local transmission infection
        
        return dSdt, dEdt, dIudt, dIcdt, dR1dt, dR2dt, dDdt, dCdt, dCum_det_dt, dCumdt, dCum_localdt
    
    # Initial conditions vector
    y0 = list(init.values())
    
    # Integrate the SIR equations over the time grid, t.
    ret = odeint(deriv, y0, t, args=(N, beta,beta2, gamma,alpha,impE,impI,pImpc,pIc,pIu_Ic,dr,cr))
    df = pd.DataFrame(ret)
    df.columns = init.keys()
    return df
    

def ordinal_seir_equations(y0,
                           t,
                           beta,
                           alpha,
                           gamma):
    """
    This is a function determining
    the SEIR system of ordinary differential equations.
    """
    
    S, E, I, R, C = y0
    dS = -beta*S*I
    dE = beta*S*I - alpha*E
    dI = alpha*E - gamma*I
    dR = gamma*I
    dC = beta*S*I
    return dS, dE, dI, dR, dC

def ordinal_seir(y0,
                 t,
                 beta,
                 alpha,
                 gamma):
    """
    SEIR ODE initial values solver
    """
    t = np.arange(0,t,1)
    arg = (beta,alpha,gamma)
    y = odeint(ordinal_seir_equations,
               y0,
               t,
               arg)
    df = pd.DataFrame(y)
    df.columns = ['S','E','I','R','C']
    return df


def seir_Rn_equations(y0,
                      t,
                      Rn,
                      dInc,
                      dInf):
    """
    Modified version of SEIR models with effective reproductive number
    """
    S, E, I, R, CI, = y0
    dS = -(1/dInf)*Rn*I
    dE = (1/dInf)*Rn*I - (1/dInc)*E
    dI = (1/dInc)*E - (1/dInf)*I
    dR = (1/dInf)*I
    dCI = (1/dInc)*E
    return dS, dE, dI, dR, dCI

def seir_Rn_solver(y0,
                   t,
                   Rn,
                   dInc,
                   dInf):
    """
    SEIQR Solver
    """
    args = (Rn,
            dInc,
            dInf)

    t = np.arange(0,t,1)
    df = pd.DataFrame(
        odeint(seir_Rn_equations,y0,t,args))

    df = df.reset_index()
    df.columns = ['DAYS','S','E','I','R','CI']
    dCI = df.loc[:,['CI']].diff()
    dCI.columns = ['dCI']
    df = pd.concat([df,dCI],axis=1)

    epicurve = df.loc[:,['dCI']]
    epicurve['week'] = np.floor(epicurve.index/7)
    epicurve = epicurve.groupby('week')['dCI'].sum()
    
    return df,epicurve
