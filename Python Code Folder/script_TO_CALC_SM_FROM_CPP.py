import os
import sys
import math
import timeit
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d, Delaunay, delaunay_plot_2d
import re
import glob
import argparse

def main(dir_name):
    file_list = glob.glob(dir_name + "*P.txt")
    for file in file_list:
        file_samp = file
        filetext= open(file_samp, 'r')
        content = filetext.read()
        filetext.close()
        lines = content.split('\n')

        normfactor_index = content.find("normfactor = \n")
        d2edgamma2_index = content.find("d2Edgamma2 = \n")
        fullboundaryvector_index = content.find("Full Boundary Vector = \n")
        hessian_index = content.find("\nHessian = \n")

        normfactor_lines = content[normfactor_index+len('\nnormfactor = \n')-1: d2edgamma2_index-1]
        d2edgamma2_lines = content[d2edgamma2_index + len('\nd2Edgamma2 = \n') -1 : fullboundaryvector_index-1]
        fullboundaryvector_lines = content[fullboundaryvector_index + len("\nFull Boundary Vector = \n") - 1 : hessian_index - 1]
        hessian_lines = content[hessian_index+len("\nHessian = "):]
        try:
            fullboundaryvector = np.float64(np.char.strip(np.array(fullboundaryvector_lines.split('\n'))))
            normfactor = np.float64(normfactor_lines)
            d2edgamma2 = np.float64(d2edgamma2_lines)

            hessian_vec = hessian_lines.split('\n')
            while ("" in hessian_vec):
                hessian_vec.remove("")

            hessian_string_array = [lines.split(' ') for lines in hessian_vec]
            for subline in hessian_string_array:
                while ("" in subline):
                    subline.remove("")
            hessian = np.float64(hessian_string_array)

            eigenvalues = np.linalg.eig(hessian).eigenvalues
            eigenvectors = np.linalg.eig(hessian).eigenvectors

            SM_affine = normfactor*d2edgamma2
            SM_nonaffine = normfactor*np.sum([1 / eigenvalues[i] * (np.dot(fullboundaryvector,eigenvectors[:,i]))**2 for i in range(len(eigenvalues)) if eigenvalues[i] > 1e-4 ])
            SM = SM_affine - SM_nonaffine

            writefile = open(file_samp[:-4] + '_PSM.txt', 'w')
            writestring = 'affine = ' + str(SM_affine) +'\n' +'nonaffine = ' + str(SM_nonaffine) + '\n' + 'SM = ' + str(SM) + '\n'
            writefile.write(writestring)
        except:
            print(file + " failed for some reason.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dir_name")
    args = parser.parse_args()
    main(args.dir_name)