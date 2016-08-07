import numpy as np
import matplotlib.pyplot as plt

def get_results():
    py_file = 'spheres_python_times.txt'
    fr_file = 'spheres_fortran_times.txt'

    py_res = {}
    fr_res = {}


    with open(py_file, 'r') as f:
        line = f.readline()
        Ncols = int(line.split()[1])
        for _ in range(3):
            line = f.readline()
            rv = float(line.split()[1])
            py_res[rv] = {}
            for _ in range(5):
                line = f.readline()
                Ns = int(line.split()[1])
                line = f.readline()
                tm = float(line.split()[2])
                py_res[rv][Ns] = tm


    with open(fr_file, 'r') as f:
        line = f.readline()
        Ncols = int(line.split()[4])
        line = f.readline()
        for _ in range(3*5):
            line = f.readline()
            Ns = int(line.split()[3])
            rv = float(line.split()[-1][: 4])
            if rv not in fr_res:
                fr_res[rv] = {}
            line = f.readline()
            tm = float(line.split()[3])
            fr_res[rv][Ns] = tm


    print py_res
    print fr_res

    colors = ['orange','g','b','r','y','k']
    i=0


    for rv in py_res.keys():
        x = np.array(py_res[rv].keys())
        a = np.argsort(x)
        y = np.array(py_res[rv].values())

        # plt.plot(np.log(x[a]), y[a], label='python rv={}'.format(rv), c=colors[i])
        plt.plot(x[a], y[a], label='python rv={}'.format(rv), c=colors[i])
        i+=1

    for rv in fr_res.keys():
        x = np.array(fr_res[rv].keys())
        a = np.argsort(x)
        y =  np.array(fr_res[rv].values())

        # plt.plot(np.log(x[a]), y[a], label='fortran rv={}'.format(np.ceil(rv * 10) / 10.), c=colors[i])
        plt.plot(x[a], y[a], label='fortran rv={}'.format(np.ceil(rv*10)/10.), c=colors[i])
        i += 1

    plt.ylabel("Runing time [secs]")
    plt.xlabel("Log Number of spheres")
    plt.legend(bbox_to_anchor=(0.5, 1))
    plt.show()




get_results()