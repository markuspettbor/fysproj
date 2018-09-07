import matplotlib.pyplot as plt

def pretty(xvar, yvar, xname, yname, header):
    # Creates single, pretty pretty plot
    plt.plot(xvar, yvar, label = yname, linewidth  = 0.8)
    plt.xlabel(xname)
    plt.ylabel(yname)
    plt.title(header)
    plt.show()

def sub_pretty(xvar, yvar1, yvar2, xname, yname1, yname2, t1, t2 = ''):
    # Creates pretty pretty subplot for two vars.
    plt.subplot(211)
    plt.title(t1)
    plt.plot(xvar, yvar1, '-k', linewidth = 0.8)
    plt.ylabel(yname1)
    plt.subplot(212)
    plt.plot(xvar, yvar2, '-c', linewidth = 0.8)
    plt.ylabel(yname2)
    plt.xlabel(xname)
    plt.title(t2)
    plt.tight_layout()
    plt.show()
