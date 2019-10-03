import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')


def boxplot(L, xticlabel, xlabel, geneName, out_file_name):
    '''
    Given a numerical array, make a box plot and save as png file
    ---
    Input: A
    An array containing numerical values
    Input: out_file_name
    Output file name
    ---
    Output
    A png file with the plot
    '''
    width = 10
    height = 3
    fig = plt.figure(figsize=(width, height), dpi=300)
    ax = fig.add_subplot(1, 1, 1)
    ax.boxplot(L)
    ax.set_xticklabels(xticlabel, rotation=90)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Gene read counts")
    ax.set_title('Gene: ' + geneName)

    plt.savefig(out_file_name, bbox_inches='tight')
    plt.close()

    return 0
