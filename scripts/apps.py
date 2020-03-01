from flask import Flask,render_template,jsonify, Response, request
from waitress import serve
import pandas as pd
import numpy as np
from models import seiqr_solver,to_beta
import matplotlib.pyplot as plt
import json

def seiqr(model,
          y0=[1000000,0,1,0,0,0,0], # [S,E,I,Q,R,CI,CQ]
          t=365,
          R0_I=2.2,
          R0_Q=0,
          epsilon_E=0.2,
          epsilon_I=0.01,
          dInc=5.2,
          dInf=2.3):
    """
    Loading model configuration
    """
    filename = f'../params/params.json'
    with open(filename,'r') as f:
        params = json.loads(f.read())
    for param in params:
        if param['model'] == model:
            selected_model = param

    """
    Parameters configuration
    """
    alpha = 1/dInc
    gamma = 1/dInf

    """
    Model solver
    """
    #y0 = [N,0,I0,0,0,0,0]
    N = np.array(y0)[:-2].sum()
    beta_I = to_beta(N=N,durations=dInf,R0=R0_I)
    beta_Q = to_beta(N=N,durations=dInf,R0=R0_Q)

    resultDays, resultWk = seiqr_solver(y0=y0,
                                        t=t,
                                        beta_I=beta_I,
                                        beta_Q=beta_Q,
                                        alpha=alpha,
                                        epsilon_E=epsilon_E,
                                        epsilon_I=epsilon_I,
                                        gamma=gamma)
    return resultDays, resultWk

app = Flask(__name__, template_folder="templates")


# Create a URL route in our application for "/"
@app.route('/')
def home():
    y0 = json.loads(request.args.get('y0','[1000000,0,1,0,0,0,0]'))
    y0 = list(map(lambda x: int(x),y0))
    return str(y0)

@app.route('/seiqr')
def seiqrapi():
    y0 = json.loads(request.args.get('y0','[1000000,0,1,0,0,0,0]'))
    t = int(request.args.get('t',365))
    R0_I = float(request.args.get('R0_I',2.2))
    R0_Q = float(request.args.get('R0_Q',0))
    epsilon_E = float(request.args.get('epsilon_E',0.01))
    epsilon_I = float(request.args.get('epsilon_I',0.2))
    dInc = float(request.args.get('dInc',5.2))
    dInf = float(request.args.get('dInf',2.3))

    model_result_days, model_result_weeks = seiqr(model='seiqr1',
                                                  y0=y0,
                                                  #N=N,
                                                  t=t,
                                                  #I0=I0,
                                                  R0_I=R0_I,
                                                  R0_Q=R0_Q,
                                                  epsilon_I=epsilon_I,
                                                  epsilon_E=epsilon_E,
                                                  dInc=dInc,
                                                  dInf=dInf)
    model_result_days = model_result_days.fillna(0)
    model_result_weeks = model_result_weeks.fillna(0)
    resultsDays = model_result_days.iloc[:,1:].values.tolist()
    resultsCols = list(model_result_days.columns)
    resultsWeeks = model_result_weeks.values.tolist()

    params = {'y0':y0,
              #'I0':I0,
              'R0_I':R0_I,
              'beta_I':to_beta(N=y0[0],durations=dInf,R0=R0_I),
              'R0_Q':R0_Q,
              'beta_Q':to_beta(N=y0[0],durations=dInf,R0=R0_Q),
              'epsilon_E':epsilon_E,
              'epsilon_I':epsilon_I,
              'dInc':dInc,
              'dInf':dInf}

    summary = {'TotalRecovered':model_result_days['R'].iloc[-1],
               'TotalInfected':model_result_days['CI'].iloc[-1],
               'TotalQuarantine':model_result_days['CQ'].iloc[-1]}
    #return Response(model_result_days.to_json(orient='records'),
    #                mimetype='application/json')
    return jsonify({'dSummary':summary,
                    'Params':params,
                    'ResultsDays':resultsDays,
                    'ResultsDaysCols':resultsCols,
                    'NewCaseWeeks':resultsWeeks})

@app.route('/seiqrplt')
def seiqrplt():
    N = int(request.args.get('N',10000000))
    t = int(request.args.get('t',365))
    I0 = int(request.args.get('I0',1))
    R0_I = float(request.args.get('R0_I',2.2))
    R0_Q = float(request.args.get('R0_Q',2.2))
    epsilon_E = float(request.args.get('epsilon_E',0))
    epsilon_I = float(request.args.get('epsilon_I',0))
    dInc = float(request.args.get('dInc',5.2))
    dInf = float(request.args.get('dInf',2.3))

    model_result_days, model_result_weeks = seiqr(model='seiqr1',
                                                  N=N,
                                                  t=t,
                                                  I0=I0,
                                                  R0_I=R0_I,
                                                  R0_Q=R0_Q,
                                                  epsilon_I=epsilon_I,
                                                  epsilon_E=epsilon_E,
                                                  dInc=dInc,
                                                  dInf=dInf)
    fig, ax = plt.subplots(1,3,figsize=(20,5))
    ax = ax.flat
    ax[0].plot(model_result_days['CI'])
    ax[1].plot(model_result_days['dCI'])
    ax[2].plot(model_result_weeks)
    fig.tight_layout()
    return fig

if __name__ == '__main__':
    #app.run(debug=True)
    serve(app)
