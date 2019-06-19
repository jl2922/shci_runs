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
    return energy, uncert

def plot_shci():
    rs = [1.50, 1.55, 1.60, 1.68, 1.80, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25]

    es = []
    e_atom = -1049.7459704502
    for r in rs:
        result_file = '2z_12e/r%.2f/OPT_eps1_2e-4_Arrow_CAS_core/result.json' % r
        energy, uncert = get_extrapolated_energy(result_file, 5)
        es.append((energy - e_atom * 2) * HA2EV)
    f = interp1d(rs, es, kind='cubic')
    x_fit = np.linspace(1.5, 3.25, num=500)
    curve_2z_12e, = plt.plot(x_fit, f(x_fit), color='orange', linestyle='dashed')
    plt.plot(rs, es, color='orange', marker='o', linestyle='none', alpha=0.5)

    es = []
    e_atom = -1049.7603535316
    for r in rs:
        result_file = '3z_12e/r%.2f/OPT_eps1_5e-4_CAS_core/result.json' % r
        energy, uncert = get_extrapolated_energy(result_file, 5)
        es.append((energy - e_atom * 2) * HA2EV)
    f = interp1d(rs, es, kind='cubic')
    x_fit = np.linspace(1.5, 3.25, num=500)
    curve_3z_12e, = plt.plot(x_fit, f(x_fit), color='green', linestyle='dashed')
    plt.plot(rs, es, color='green', marker='o', linestyle='none', alpha=0.5)

    es = []
    e_atom = -1049.9325698806
    for r in rs:
        result_file = '2z_28e/r%.2f/result.json' % r
        energy, uncert = get_extrapolated_energy(result_file, 6)
        es.append((energy - e_atom * 2) * HA2EV)
    f = interp1d(rs, es, kind='cubic')
    x_fit = np.linspace(1.5, 3.25, num=500)
    curve_2z_28e, = plt.plot(x_fit, f(x_fit), color='blue', linestyle='dashed')
    plt.plot(rs, es, color='blue', marker='o', linestyle='none', alpha=0.5)

    es = []
    e_atom = -1049.9325698806
    for r in rs:
        result_file = '2z_28e_HFC/r%.2f/result.json' % r
        energy, uncert = get_extrapolated_energy(result_file, 6)
        es.append((energy - e_atom * 2) * HA2EV)
    f = interp1d(rs, es, kind='cubic')
    x_fit = np.linspace(1.5, 3.25, num=500)
    curve_2z_28e_hfc, = plt.plot(x_fit, f(x_fit), color='blueviolet', linestyle='dashed')
    plt.plot(rs, es, color='blueviolet', marker='o', linestyle='none', alpha=0.5)

    return curve_2z_12e, curve_3z_12e, curve_2z_28e, curve_2z_28e_hfc

def plot():
    plt.figure(figsize=(5.5, 4.0))
    params = {'mathtext.default': 'regular' }
    plt.rcParams.update(params)
    curve_experiment = plot_experiment()
    curve_2z_12e, curve_3z_12e, curve_2z_28e, curve_2z_28e_hfc = plot_shci()
    #curve_uhf, curve_ccsd, curve_ccsdt = plot_hfcc()

    plt.xlabel('Bond Length ($\AA$)')
    plt.ylabel('Atomization Energy (eV)')
    plt.title('Cr$_2$ Potential Energy Curve')
    curves = (curve_2z_12e, curve_3z_12e, curve_2z_28e, curve_2z_28e_hfc, curve_experiment)
    labels = ('SHCI 12e cc-pVDZ (CAS Core)', 'SHCI 12e cc-pVTZ (CAS Core)', 'SHCI 28e cc-pVDZ (CAS Core)', 'SHCI 28e cc-pVDZ (HF Core)', 'Experiment')
    plt.legend(curves, labels)
    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)
    plt.tight_layout()
    plt.grid(True, ls=':')
    plt.savefig('cr2curve.png', format='png', dpi=300)
    plt.show()


if __name__ == '__main__':
    plot()
