import matplotlib
#matplotlib.use('TkAgg')
from matplotlib import pyplot as plt

import seaborn as sns
from statistics import mean

#sns.set()

def tr_makeplot(trace, xlabel="Sample [Pt]", ylabel="Power Consumption [mV]"):
    fig = plt.figure()
    if type(trace[0]) is int:
        coeff = 1
    else:
        coeff = 1000
    plt.plot([coeff*s for s in trace])
    plt.grid(True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()

def avg_tr_makeplot(tab, xlabel="Sample [Pt]", ylabel="Power Consumption [mV]", title="Average Power Consumption"):
    fig = plt.figure()
    if type(tab[0][0]) is int:
        coeff = 1
    else:
        coeff = 1000
    plt.plot([coeff*mean(l) for l in map(list, zip(*tab))])
    plt.grid(True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()

def avg_fft_makeplot(freq, tabf, xlabel="Frequency [Hz]", ylabel="Fourier Transform [mV s]", title="Average Power Consumption\nunder Fourier Transform"):
    fig = plt.figure()
    if type(tabf[0][0]) is int:
        coeff = 1
    else:
        coeff = 1000
    plt.plot(freq, [coeff*mean(l) for l in map(list, zip(*tabf))])
    plt.grid(True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()

def cmp_avg_makeplot(first, second, xlabel="Sample [Pt]", ylabel="Power Consumption [mV]", title="Average Power Consumption"):
    fig = plt.figure()
    plt.title(title)
    if type(first[0][0]) is int:
        coeff = 1
    else:
        coeff = 1000
    plt.subplot(2, 1, 1)
    plt.plot([coeff*mean(l) for l in map(list, zip(*first))])
    plt.grid(True)
    plt.xlabel("Sample [Pt]")
    plt.ylabel("mV")
    plt.subplot(2, 1, 2)
    if type(second[0][0]) is int:
        coeff = 1
    else:
        coeff = 1000
    plt.plot([coeff*mean(l) for l in map(list, zip(*second))])
    plt.grid(True)
    plt.xlabel("Sample [Pt]")
    plt.ylabel("mV")
    plt.show()

def pcc_vs_pcc_makeplot(first, second, xlabel="Sample [Pt]", ylabel="Power Consumption [mV]", title="Correlation Power Analysis Comparison"):
    fig = plt.figure()
    plt.title(title)
    if type(first[0][0]) is int:
        coeff = 1
    else:
        coeff = 1000
    plt.subplot(2, 1, 1)
    plt.plot([coeff*mean(l) for l in map(list, zip(*first))])
    plt.grid(True)
    plt.xlabel("Sample [Pt]")
    plt.ylabel("mV")
    plt.subplot(2, 1, 2)
    if type(second[0][0]) is int:
        coeff = 1
    else:
        coeff = 1000
    plt.plot([coeff*mean(l) for l in map(list, zip(*second))])
    plt.grid(True)
    plt.xlabel("Sample [Pt]")
    plt.ylabel("mV")
    plt.show()

def pcc_makeplot(pcc, xlabel="Sample [Pt]", ylabel="PCC", title="Correlation Power Analysis", fit=True):
    fig = plt.figure()
    plt.plot(pcc)
    plt.grid(True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if not fit:
        ytop = max(max(pcc), 0.5)
        ybot = min(min(pcc), -0.5)
        plt.ylim(ybot, ytop)
    plt.title(title)
    plt.show()

def pcc_vs_avg_tr_makeplot(tr, pcc, title="Correlation Location"):
    fig = plt.figure()
    plt.subplot(2, 1, 1)
    plt.title(title)
    plt.plot([1000*mean(l) for l in map(list, zip(*tr))])
    plt.grid(True)
    plt.xlabel("Sample [Pt]")
    plt.ylabel("mV")
    plt.subplot(2, 1, 2)
    plt.plot(pcc, color='green')
    plt.grid(True)
    plt.xlabel("Sample [Pt]")
    plt.ylabel("PCC")
    plt.show()

def pcc_hyp_makeplot(tab, xlabel="Sample [Pt]", ylabel="PCC", title="Correlation Power Analysis", fit=True):
    fig = plt.figure()
    for i in range(len(tab)):
        plt.plot(tab[i], label=f'hyp. #{i+1}')
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.grid(True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if not fit:
        ytop = max(max(max(tab)), 0.5)
        ybot = min(min(min(tab)), -0.5)
        plt.ylim(ybot, ytop)
    plt.title(title)
    plt.tight_layout()
    plt.show()

def pcc_2hyp_makeplot(tab, xlabel="Sample [Pt]", ylabel="PCC", title="Correlation Power Analysis", fit=True):
    fig = plt.figure()
    plt.plot(tab[0], label=f'valid', color='green')
    plt.plot(tab[1], label=f'wrong', color='red')
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.grid(True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if not fit:
        ytop = max(max(max(tab)), 0.5)
        ybot = min(min(min(tab)), -0.5)
        plt.ylim(ybot, ytop)
    plt.title(title)
    plt.tight_layout()
    plt.show()

def fft_pcc_hyp_makeplot(freq, tabf, xlabel="Frequency [Hz]", ylabel="PCC", title="Correlation Power Analysis\nunder Fourier Transform", fit=True):
    fig = plt.figure()
    for i in range(len(tabf)):
        plt.plot(freq, tabf[i], label=f'hyp. #{i+1}')
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.grid(True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if not fit:
        ytop = max(max(max(tabf)), 0.5)
        ybot = min(min(min(tabf)), -0.5)
        plt.ylim(ybot, ytop)
    plt.title(title)
    plt.tight_layout()
    plt.show()

def fft_pcc_2hyp_makeplot(freq, tabf, xlabel="Frequency [Hz]", ylabel="PCC", title="Correlation Power Analysis\nunder Fourier Transform", fit=True):
    fig = plt.figure()
    plt.plot(freq, tabf[0], label=f'valid', color='green')
    plt.plot(freq, tabf[1], label=f'wrong', color='red')
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.grid(True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if not fit:
        ytop = max(max(max(tabf)), 0.5)
        ybot = min(min(min(tabf)), -0.5)
        plt.ylim(ybot, ytop)
    plt.title(title)
    plt.tight_layout()
    plt.show()

def multi_traces_makeplot(tab, xlabel="Sample [Pt]", ylabel="Power Consumption [mV]", title="Power Consumption"):
    fig = plt.figure()
    if type(tab[0][0]) is int:
        coeff = 1
    else:
        coeff = 1000
    for i in range(len(tab)):
        plt.plot([coeff*t for t in tab[i]], label=f'trace #{i+1}')
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.grid(True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.show()

def multi_traces_zoom_makeplot(tab, xlabel="Sample [Pt]", ylabel="Power Consumption [mV]", title="Power Consumption"):
    fig = plt.figure()
    plt.subplot(1, 2, 1)
    if type(tab[0][0]) is int:
        coeff = 1
    else:
        coeff = 1000
    for i in range(len(tab)):
        plt.plot([coeff*t for t in tab[i]], label=f'trace #{i+1}')
    plt.grid(True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.subplot(1, 2, 2)
    for i in range(len(tab)):
        plt.plot(np.linspace(190, 219, num=30), [coeff*t for t in tab[i][190:220]], label=f'trace #{i+1}')
    plt.grid(True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()