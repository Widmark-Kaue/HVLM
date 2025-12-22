import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def set_aiaa_style(): # type: ignore
    """
    AIAA-style matplotlib preset.
    
    single_column = True  -> figura de 1 coluna
    single_column = False -> figura de 2 colunas
    """
    # if figsize == None:
    #     if single_column:
    #         figsize = (3.5, 2.5)   # ~8.9 cm (1 coluna)
    #     else:
    #         figsize = (7.0, 4.5)   # ~17.8 cm (2 colunas)

    plt.rcParams.update({

        # # === Figura ===
        # 'figure.dpi': 300,
        # 'savefig.dpi': 600,
        # 'savefig.bbox': 'tight',
        # 'savefig.pad_inches': 0.02,

        # === Fonte ===
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif'],
        'font.size': 8,

        # === Eixos ===
        'axes.labelsize': 8,
        'axes.titlesize': 9,
        'axes.linewidth': 0.8,
        'axes.grid': True,
        'axes.axisbelow': True,

        # === Grid ===
        'grid.linestyle': ':',
        'grid.linewidth': 0.5,
        'grid.alpha': 0.5,

        # === Linhas ===
        'lines.linewidth': 1.5,
        'lines.markersize': 4,
        'lines.markeredgewidth': 0.8,

        # === Ticks ===
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.major.size': 4,
        'ytick.major.size': 4,
        'xtick.minor.size': 2,
        'ytick.minor.size': 2,

        # === Legenda ===
        'legend.fontsize': 7,
        'legend.frameon': False,
        'legend.handlelength': 2.0,

        # === MathText ===
        'mathtext.fontset': 'cm',
        'mathtext.rm': 'serif',

    })

def plot_polar(vlmData:np.ndarray, path_case:Path, savefig:bool = False):
    #paths
    path_images = Path('images')
    case = path_case.parent.name
    
    # data
    alpha = vlmData[:, 0]
    CL = vlmData[:, 1]
    CDi = vlmData[:, 2]
    avlData = np.loadtxt(path_case)

    plt.figure(figsize=(8, 3.5))
    plt.subplot(1, 2, 1)
    plt.plot(avlData[:,0], avlData[:,1], 'k-o', label = 'AVL')
    plt.plot(alpha, CL, 'r--s', label = 'code')
    # plt.yticks(avlData[:,1])
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{L}$')
    # plt.grid()
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(avlData[:,0], avlData[:,2], 'k-o', label = 'AVL')
    plt.plot(alpha, CDi, 'r--s', label = 'code')
    # plt.plot(alpha, CDi2, 'm--^', label = 'code 2')
    # plt.yticks(avlData[:,1])
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{D_i}$')
    # plt.legend()

    if savefig:
        plt.tight_layout()
        plt.savefig(path_images.joinpath(f'{case}.pdf'), dpi = 600, format = 'pdf')
    plt.suptitle(case.upper())
    plt.tight_layout()
    plt.show(block = False)