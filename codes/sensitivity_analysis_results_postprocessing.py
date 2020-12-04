#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import h5py as hp

params = {'1':'a',
          '2':'D',
          '3':'A',
          '4':'R',
          '5':'\chi',
          '6':'A_{\sigma}',
          '7':'K',
          '8':'\gamma'
          }

# interval = 'IC3'
# indices = open('Sensitivity_data/PCE_data/data_{0}/Pygpc_Sobol_idx.txt'.format(interval), 'r')
filepath = 'Sensitivity_data/'
indices = open(filepath+'Pygpc_Sobol_idx.txt', 'r')

# indices = open('Sensitivity_data/Pygpc_Sobol_idx.txt', 'r')
first_order = []
scnd_order = []
total_order = []
while True:

    # Get next line from file
    line = indices.readline()
    tmp = line.split(',')
    # print(tmp)
    if len(tmp[0])==1:
        first_order.append([tmp[0], np.float64(tmp[1])])
    if len(tmp[0])==3:
        tmp_1 = tmp[0].split(' ')
        scnd_order.append([tmp_1[0],tmp_1[1], np.float64(tmp[1])])
    # if line is empty
    # end of file is reached
    if not line:
        break
first_order = np.asarray(first_order)
scnd_order = np.asarray(scnd_order)
indices.close()

def get_total_idx(par):
    idx = 0
    for i in range(len(scnd_order[:,0:2])):
        if np.array2string(scnd_order[:,0:2][i]).find(par)>-1:
            idx = idx + np.float64(scnd_order[:,2][i])

    f_idx = np.int(np.where(first_order[:, 0] == par)[0])
    idx = idx + np.float64(first_order[f_idx, 1])
    return idx

for id in first_order[:,0]:
    total_order.append([id, get_total_idx(id)])

total_order = np.asarray(total_order)
# indices = np.loadtxt('Sensitivity_data/Pygpc_Sobol_idx.txt', delimiter=',', converters={0:str, 1: np.float64}, dtype=(str, np.float64))

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['font.sans-serif'] = ['Tahoma']
plt.rcParams['figure.figsize'] = [8., 8.]
plt.rcParams.update({'font.size': 18})


params_list_plt = []
for i in range(len(first_order)):
    params_list_plt.append('$'+params.get(first_order[i,0])+'$')

pos = np.arange(len(first_order))  # the label locations
width = 0.35  # the width of the bars

color1 = '#1f77b4'
color2 = '#ff7f0e'
# color2 = '#ac1129'
color2 = '#ececec'


#  --------------------------------- Plot first and total order indices ---------------------------
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
rects1 = ax1.bar(pos - width/2,
                 np.abs(np.float64(first_order[:,1])),
                 capsize=5,
                 width=width,
                 fill= False,
                 edgecolor='black',
                 linestyle = '-',
                 linewidth=3,
                 label='1st Sob. ind.')
rects2 = ax2.bar(pos + width/2,
                 np.abs(np.float64(total_order[:,1])),
                 capsize=5, width=width,
                 fill= False, hatch='//\\',
                 color='k',
                 edgecolor='black',
                 linewidth=3,
                 label='Tot. Sob. ind.')

ax1.set_ylabel('First order Sobol indices', color='k')
ax2.set_ylabel('Total order Sobol indices', color='black')
# ax.set_title('Sobol indices for $\mu_1$')
ax1.set_ylim(ymin=0.,ymax=np.max(np.abs(np.float64(first_order[:,1])))+np.min(np.abs(np.float64(first_order[:,1]))))
ax1.set_xticks(pos)
ax1.set_xticklabels(params_list_plt)
ax1.tick_params(axis='y', labelcolor='k')
ax2.tick_params(axis='y', labelcolor='black')
ax1.legend()

ax1.legend(loc='upper left')
ax2.legend(loc='upper right')

#  --------------------------------- Plot second order indices ---------------------------
# Second order param lists
scnd_params_list = []
for j in range(len(scnd_order)):
    scnd_params_list.append('${0}$, ${1}$'.format(params.get(scnd_order[:,0:2][j][0]), params.get(scnd_order[:,0:2][j][1])))

fig = plt.figure(2)
size_of_fig=6
# fig.set_size_inches(size_of_fig+2, size_of_fig)
# plt.bar(np.arange(len(scnd_params_list)), scnd_indices, width=.8, yerr=err2, capsize=7, edgecolor='black', linewidth=3);
plt.bar(np.arange(len(scnd_order)),
        np.abs(np.float64(scnd_order[:,2])),
        width=.8,
        fill= False,
        hatch='//\\',
        color='k',
        edgecolor='black',
        capsize=7,
        linewidth=3);
plt.ylabel('Second order Sobol indices');
# plt.title('Sobol second order indices for $\mu_1$');
# plt.ylim(ymin=0.,ymax=0.3)
plt.xticks(np.arange(len(scnd_order)), scnd_params_list)
plt.xticks(rotation=60)



# ---------------------------------------------------------------------------------------------------------------------
# -------------------------------- Validate the gpc approximation against the original model --------------------------
# ---------------------------------------------------------------------------------------------------------------------

# validation_plot = 'Sensitivity_data/PCE_data/data_{0}/PCE_dataplot.hdf5'.format(interval)
validation_plot = filepath+'PCE_dataplot.hdf5'
# validation_mc = 'Sensitivity_data/PCE_data/data_{0}/PCE_datamc.hdf5'.format(interval)
validation_mc = filepath+'PCE_datamc.hdf5'

validation_plot_data = hp.File(validation_plot, "r+")
validation_mc_data = hp.File(validation_mc, "r+")

# keys for validation_plot_data : 'gpc_vs_original_plot' >> 'grid'              >> ('coords', 'coords_norm')
#                                                        >> 'model_evaluations' >> ('difference', 'gpc', 'original', 'original_all_qoi')


size = 10
mu1_a_A_gpc = np.reshape(validation_plot_data['gpc_vs_original_plot']['model_evaluations']['gpc'], (size, size)).T
mu1_a_A_orig = np.reshape(validation_plot_data['gpc_vs_original_plot']['model_evaluations']['original'], (size, size)).T

a_s = np.reshape(validation_plot_data['gpc_vs_original_plot']['grid']['coords'][:,0], (size,size))
A_s = np.reshape(validation_plot_data['gpc_vs_original_plot']['grid']['coords'][:, 2], (size, size))

plt.figure(3)
d_delta = np.max(A_s.T[1:size, 0]) / np.float64(size)

d_mu = np.min(np.exp(mu1_a_A_gpc[0, 1:size - 2])) - np.min(np.exp(mu1_a_A_gpc[0, 1:size - 1]))
plt.plot(A_s.T[1:size, 0], np.exp(mu1_a_A_orig[0, 1:size]), color='k', linewidth=3, label='original model')
plt.plot(A_s.T[1:size, 0], np.exp(mu1_a_A_gpc[0, 1:size]), '--', color='k', linewidth=3, label='gPC approximation')


d_mu = np.min(np.exp(mu1_a_A_gpc[size-1, 1:size - 2])) - np.min(np.exp(mu1_a_A_gpc[size - 1, 1:size - 1]))
plt.plot(A_s.T[1:size, 0], np.exp(mu1_a_A_orig[size - 1, 1:size]), color='k', linewidth=3)
plt.plot(A_s.T[1:size, 0], np.exp(mu1_a_A_gpc[size - 1, 1:size]), '--', color='k', linewidth=3)
plt.xlabel('Im. strength $A$ in $cell_c^{-1} \cdot day^{-1}$')
plt.ylabel('tumor diameter in $mm$')
plt.legend()


# keys for validation_mc_data : 'gpc_vs_original_mc' >> 'error'             >> 'qoi_0'                     >> 'nrmsd'
#                                                    >> 'grid'              >>  ('coords', 'coords_norm')
#                                                    >> 'model_evaluations' >> 'results'                   >> ('gpc', 'orig')
#                                                    >> 'pdf'               >> 'qoi_0'                     >> ('gpc', 'original')

# -------------------------------- Validate the gpc pdf approximation against the original model ----------------------

# res = np.loadtxt('Sensitivity_data/back_up/results.txt')
gpc_pdf = validation_mc_data['gpc_vs_original_mc']['pdf']['qoi_0']['gpc']
orig_pdf = validation_mc_data['gpc_vs_original_mc']['pdf']['qoi_0']['original']

res_gpc = validation_mc_data['gpc_vs_original_mc']['model_evaluations']['results']['gpc']
res_orig = validation_mc_data['gpc_vs_original_mc']['model_evaluations']['results']['orig']

sns.distplot(np.exp(res_gpc), hist=True, kde=True,
             bins=int(500/5), color = 'gray',
             hist_kws={'edgecolor':'gray', 'range':(0.,0.0009)},
             kde_kws={'shade': True, 'label':'original model', 'color':'k','linestyle':'--', 'linewidth': 3,'clip':(0., 0.0009)})
sns.distplot(np.exp(res_orig), hist=True, kde=True,
             bins=int(500/5), color = 'k',
             hist_kws={"histtype": "step", 'edgecolor':'black', 'range':(0.,0.0009)},
             kde_kws={'shade': True,  'label':'gPC approximation', 'linewidth': 1, 'clip':(0., 0.0009)})
plt.xlabel('$\\mu_1$', fontsize=12)
plt.ylabel('$P(\\mu_1)$', fontsize=16)
plt.legend()

fig = plt.figure(2)
size_of_fig=8
fig.set_size_inches(size_of_fig+6, size_of_fig+2)
sns.distplot(res_orig, hist=True, kde=True,
             bins=int(100/5), color = 'gray',
             hist_kws={'edgecolor':'black', 'range':(np.min(res_gpc),np.max(res_gpc))},
             kde_kws={'shade': True, 'label':'original model', 'color':'k','linestyle':'dotted', 'linewidth': 10, 'clip':(np.min(res_gpc),np.max(res_gpc))})
sns.distplot(res_gpc,  hist=True, kde=True,
             bins=int(180/5), color = 'k',
             hist_kws={"histtype": "step", 'edgecolor':'black', 'range':(np.min(res_orig),np.max(res_orig))},
             kde_kws={'shade': True,  'label':'gPC approximation', 'linewidth': 2, 'clip':(np.min(res_orig),np.max(res_orig))})
plt.xlabel('$\\ln(\\mu_1)$', fontsize=18)
plt.ylabel('$P(\\ln(\\mu_1))$', fontsize=18)
plt.legend()