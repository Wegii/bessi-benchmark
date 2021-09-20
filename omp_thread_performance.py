import matplotlib.pyplot as pyplot
import csv
import statistics
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
from matplotlib import animation
from scipy import interpolate

class core:
    def __init__(self):
        self.version = "0.0.1"

    def get_content(files):
        contents = []

        for file in files:
            f = csv.reader(open(file), delimiter=',')

            x = []
            for row in f:
                x.append(float(row[0]))
        
            contents.append(x)

        return contents

    def make_graph(x, y, labels, title, label_x, label_y, location, legend='upper right', log=False,
                                sub_title = ''):
        fig, axes = pyplot.subplots(1, figsize=(10,4.5))

        for i in range(0, len(labels)):
            axes.plot(x[i], y[i], label=labels[i], marker = 'o', markersize=3)
            axes.set_title(title)
            axes.set_xlabel(label_x)
            axes.set_ylabel(label_y)
            if log:
                axes.set_xscale('log', base=2)

        pyplot.legend(loc=legend)
        pyplot.savefig(location, dpi=300, bbox_inches='tight')
        pyplot.show()


    def make_graph_error(x, y, error, labels, title, label_x, label_y, location, legend='upper left',
                    sub_title = ''): 
        fig, axes = pyplot.subplots(1, figsize=(10,4.5))
        
        for i in range(0, len(labels)):
            axes.plot(x[i], y[i], label=labels[i], markersize=3)
            axes.set_title(title)
            axes.set_xlabel(label_x)
            axes.set_ylabel(label_y)
            pyplot.fill_between(x[i], y[i]-error[i], y[i]+error[i], edgecolor='none', alpha=0.2)

        axes.set_xscale('log')
        fig.tight_layout()
        pyplot.margins(0, 0)
        pyplot.legend(loc=legend)
        pyplot.savefig(location, dpi=300, bbox_inches='tight')
        pyplot.show()


    def make_graph_3d(x, y, z, title, label_x, label_y, label_z, location):
        X = np.array(x).flatten()
        Y = np.array(y).flatten()
        Z = np.array(z).flatten()

        fig = pyplot.figure(figsize=(4.5,4.5))
        ax = fig.add_subplot(111, projection='3d')
                        
        ax.plot_trisurf(X, Y, Z, cmap=pyplot.cm.viridis, linewidth=0.2)

        ax.set_title(title)
        ax.set_xlabel(label_x)
        ax.set_ylabel(label_y)
        ax.set_zlabel(label_z)
        ax.view_init(elev=30, azim=120)
        pyplot.savefig(location, dpi=300)
        #pyplot.show()

    def make_heatmap(values, labels_x, labels_y, title, location):
        fig, ax = pyplot.subplots(1, figsize=(10, 7.5))
        im = ax.imshow(values)
        
        ax.set_xticks(np.arange(len(labels_x)))
        ax.set_yticks(np.arange(len(labels_y)))
        ax.set_xticklabels(labels_x)
        ax.set_yticklabels(labels_y)

        ax.set_xlabel('Chunk Size')
        ax.set_ylabel('Threads')

        pyplot.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        for i in range(len(labels_y)):
            for j in range(len(labels_x)):
                text = ax.text(j, i, ('%.2f'%values[i, j]) + '' , ha="center", va="center", color="w", fontsize=8)

        for edge, spine in ax.spines.items():
            spine.set_visible(False)

        #cbar = ax.figure.colorbar(im, ax=ax)
        cax = fig.add_axes([ax.get_position().x1 + 0.01, ax.get_position().y0, \
                                        0.02, ax.get_position().height])
        pyplot.colorbar(im, cax=cax)

        ax.set_title(title)
        pyplot.savefig(location, dpi=300, bbox_inches='tight')
       
        pyplot.show()

Core = core()
basedir = '/work/pwegmann/BESSI_TEST/bessi_test/benchmark/data/'
filedir = '/work/pwegmann/BESSI_TEST/bessi_test/benchmark/figures/'

# Execution time vs. year for different number of threads averaged 
# over the static schedule chunk size
threads = [8, 12, 24, 32, 48]
schedule = [1, 2, 4, 8, 12, 24, 48, 96]

labels = []
runtime_thread = []
x = []
error = []
for t in threads:
    runtime_schedule = []
    for s in schedule:
        content = core.get_content([basedir + "static_t" + str(t) + "_s" + str(s)])
        runtime_schedule.append(content[0])

    # Average execution time over different chunck sizes
    runtime_thread.append(np.average(runtime_schedule, axis=0))
    x.append(list(range(0, len(content[0]))))
    labels.append("T=" + str(t))

    # Calculate error
    se = []
    for j in range(0, len(runtime_schedule[0])):
        sv = []
        for s in range(0, len(schedule)):
            sv.append(runtime_schedule[s][j])
        # Get maximum deviation from average value
        se.append(max(abs(runtime_thread[-1][j] - sv))/2.0)
    
    error.append(se)

#core.make_graph_error(x, runtime_thread, error, labels, "Average execution time for different runtime configurations",
#        "Year", "Runtime [s]", filedir + "omp_thread_performance_thread_schedule_average.png")

threads = [24, 48]
schedule = [1, 2, 4, 8, 12, 24, 48]
x = []
y = []
labels = []
for t in threads:
    for s in schedule:
        content = core.get_content([basedir + "static_t" + str(t) + "_s" + str(s)])
        y.append(content[0])
        x.append(list(range(0, len(content[0]))))
        labels.append('T=' + str(t) + ", C=" + str(s))

core.make_graph(x, y, labels, "Execution time for different runtime configurations",
                "Year", "Runtime [s]", filedir + "omp_thread_performance_thread_schedule_fastest.png")


# Speedup for different number of threads and chunk size
threads = [1, 8, 12, 24, 32, 48]
schedule = [1, 2, 4, 8, 12, 24, 48, 96]
ly = ['T=' + str(i) for i in threads]
lx = ['C=' + str(i) for i in schedule]
z = []
for t in threads:
    iz = []
    for s in schedule:
        content = core.get_content([basedir + "static_t" + str(t) + "_s" + str(s)])
        # Average runtime after 80 years
        iz.append(np.average(content[0][79:]))

    z.append(iz)

# Get maximum runtime
max_z = -1
for i in range(0, len(z)):
    if i == 0:
        max_z = max(z[i])
    else:
        nm = max(z[i])
        
        if (nm > max_z):
            max_z = nm

# Speedup over maximum execution time
for i in range(0, len(z)):
    for j in range(0, len(z[0])):
        z[i][j] = max_z/z[i][j]

#core.make_heatmap(np.array(z), lx, ly, "Speedup of different runtime configurations", 
#        filedir + "omp_thread_performance_thread_schedule_heatmap.png")












