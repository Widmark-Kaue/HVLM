import matplotlib.pyplot as plt

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