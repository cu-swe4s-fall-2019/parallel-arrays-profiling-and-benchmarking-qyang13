import sys
import matplotlib
import matplotlib.pyplot as plt
import math_lib as ml
matplotlib.use('Agg')


def boxplot(A, out_file_name):
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
    plt.boxplot(A)
    plt.ylabel('Value Distribution')
    plt.xlabel('Catagory')
    plt.title('Mean: ' + str(ml.list_mean(A))
              + ' '
              + 'Stdev: ' + str(ml.list_stdev(A)))
    plt.savefig(out_file_name, bbox_inches='tight')
    pass


def histogram(A, out_file_name):
    '''
    Given a numerical array, make a histogram and save as png file
    ---
    Input: A
    An array containing numerical values
    Input: out_file_name
    Output file name
    ---
    Output
    A png file with the plot
    '''
    plt.hist(A, bins=50)
    plt.ylabel('Frequency')
    plt.xlabel('Value')
    plt.title('Mean: ' + str(ml.list_mean(A))
              + ' '
              + 'Stdev: ' + str(ml.list_stdev(A)))
    plt.savefig(out_file_name, bbox_inches='tight')
    pass


def combo(A, out_file_name):
    '''
    Given a numerical array, make a box plot and
    a histogram and save as png file
    ---
    Input: A
    An array containing numerical values
    Input: out_file_name
    Output file name
    ---
    Output
    A png file with the plots
    '''
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.boxplot(A)
    ax1.set_xlabel('Catagory')
    ax1.set_ylabel('Value Distribution')
    ax2.hist(A)
    ax2.set_xlabel('Value')
    ax2.set_ylabel('Frequency')
    fig.suptitle('Mean: ' + str(ml.list_mean(A))
                 + ' '
                 + 'Stdev: ' + str(ml.list_stdev(A)))
    fig.savefig(out_file_name, bbox_inches='tight')
    pass
