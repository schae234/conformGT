#!/usr/bin/env python3
import sys

import matplotlib
matplotlib.use('Agg')

import numpy as np 
import ponytools as pc
import matplotlib.pyplot as plt

from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram

from optparse import OptionParser


def main(args):
    vcf = pc.VCF(args.vcf).genotypes

    # Grab labels for matrix
    sample_labels = vcf.columns.values
    # cluster to samples
    sample_distance = squareform(pdist(vcf.T.as_matrix(),'correlation')) 

    # Get a figure
    plt.figure(figsize=(10,100))

    dendro = dendrogram(linkage(sample_distance,method='complete'),orientation='right',labels=sample_labels)
   
    plt.tight_layout() 
    plt.savefig(args.out)


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('--vcf', type=str)
    parser.add_option('--out', type=str)
    args,options = parser.parse_args()
    sys.exit(main(args))
