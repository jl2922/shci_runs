#!/usr/bin/env python
import argparse
import json
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import statsmodels.api as sm
import scipy.stats

parser = argparse.ArgumentParser()
parser.add_argument('--result_file', default='result.json')
parser.add_argument('--save_figure', type=bool, default=True)
parser.add_argument('--show_figure', type=bool, default=True)
parser.add_argument('--order', type=int, default=2)
args = parser.parse_args()

HA2EV = 27.2114

def plot_experiment():
    df = pd.read_csv('experiment.csv')
    r = df['r'].values
    e = df['e'].values
    f = interp1d(r, e, kind='cubic')
    x_fit = np.linspace(1.5, 3.25, num=500)
    curve, = plt.plot(x_fit, f(x_fit), color='grey', linestyle='dashed')
    return curve

def model_aug(x):
    x_aug = (x, )
    for i in range(2, (args.order + 1)):
        x_aug = x_aug + (x**i, )
    x_aug = np.column_stack(x_aug)
    x_aug = sm.add_constant(x_aug)
    return x_aug

def get_extrapolated_energy(result_file, n_points):
    result = open(result_file).read()
    result = json.loads(result)
    x = []
    y = []
    e = []
    energy_vars = result['energy_var']
    for eps_var, energy_var in energy_vars.items():
        if result['energy_total'].get(eps_var) is None:
            continue
        energy_totals = result['energy_total'][eps_var]
        eps_pt = 1.0
        for eps_pt_iter, energy_total_iter in energy_totals.items():
            eps_pt_iter = float(eps_pt_iter)
            if eps_pt_iter < eps_pt:
                eps_pt = eps_pt_iter
                energy_total = energy_total_iter['value']
        y.append(energy_total)
        x.append(energy_var - energy_total)
        e.append(float(eps_var))
    x = np.array(x)
    y = np.array(y)
    e = np.array(e)
    smallest_points = e.argsort()[:n_points]
    x = x[smallest_points]
    y = y[smallest_points]
    x_aug = model_aug(x)
    weights = 1.0 / x**2
    alpha = 0.05
    point = np.zeros(x_aug.shape[1])
    point[0] = 1.0
    t = scipy.stats.t.ppf((2 - alpha) / 2., x.shape[0] - 3)
    tt = t * 2
    fit = sm.WLS(y, x_aug, weights).fit()
    predict = fit.get_prediction(point).summary_frame(alpha=alpha)
    predict = predict.iloc[0]
    energy = fit.params[0]
    uncert = predict['mean_ci_upper'] - predict['mean']
    if np.isnan(uncert):
        uncert = 9999
    #print('(%.2f Conf.) Extrapolated Energy: %.10f +- %.10f' % ((1.0 - alpha, energy, uncert)))
    return energy, uncert, y

def plot_raw(rs, ys, index, color):
    yis = []
    for i in range(len(rs)):
        yis.append(ys[i][index])
    f = interp1d(rs, yis, kind='cubic')
    x_fit = np.linspace(1.5, 3.25, num=500)
    curve, = plt.plot(x_fit, f(x_fit), color=color, linestyle='solid')
    plt.plot(rs, yis, color=color, marker='o', linestyle='none', alpha=0.7)
    return curve

def plot_shci():
    rs = [1.50, 1.55, 1.60, 1.68, 1.80, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25]

    es = []
    e_atom = -1049.9325698806
    ys = []
    for r in rs:
        result_file = '2z_28e_HFC/r%.2f/result.json' % r
        energy, uncert, y = get_extrapolated_energy(result_file, 6)
        es.append((energy - e_atom * 2) * HA2EV)
        ys.append((y - e_atom * 2) * HA2EV)
    f = interp1d(rs, es, kind='cubic')
    x_fit = np.linspace(1.5, 3.25, num=500)
    curve_2z_28e, = plt.plot(x_fit, f(x_fit), color='blue', linestyle='solid')
    plt.plot(rs, es, color='blue', marker='o', linestyle='none', alpha=0.7)

    curve_s0 = plot_raw(rs, ys, 0, 'green')
    curve_s1 = plot_raw(rs, ys, 1, 'green')
    curve_s2 = plot_raw(rs, ys, 2, 'green')

    return curve_2z_28e, curve_s0, curve_s1, curve_s2

def plot():
    plt.figure(figsize=(5.5, 4.0))
    params = {'mathtext.default': 'regular' }
    plt.rcParams.update(params)
    curve_experiment = plot_experiment()
    curve_2z_28e, curve_s0, curve_s1, curve_s2 = plot_shci()

    plt.xlabel('Bond Length ($\AA$)')
    plt.ylabel('Atomization Energy (eV)')
    plt.title('Cr$_2$ Potential Energy Raw Curves with 28e and cc-pVDZ')
    curves = (curve_s0, curve_s1, curve_s2, curve_2z_28e, curve_experiment)
       # 'SHCI $\epsilon_1$=2e-5',
       # 'SHCI $\epsilon_1$=1e-5',
       # 'SHCI $\epsilon_1$=5e-6',
    labels = (
        'SHCI $\epsilon_1=2 \\times 10^{-5}$',
        'SHCI $\epsilon_1=1 \\times 10^{-5}$',
        'SHCI $\epsilon_1=5 \\times 10^{-6}$',
        'SHCI Extrapolated',
        'Experiment')
    plt.legend(curves, labels)
    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)
    plt.tight_layout()
    plt.grid(True, ls=':')
    plt.savefig('cr2raw.png', format='png', dpi=300)
    plt.show()


if __name__ == '__main__':
    plot()
