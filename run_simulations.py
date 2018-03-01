#!/usr/bin/env python

from subprocess import check_call
from Bio import SeqIO
import os
from joblib import Parallel, delayed
import multiprocessing

true_tree = 'sim.tree'
simultron = './physher/bin/simultron'
beast_dir = '/Applications/BEAST 2.4.5/bin'

data = {
    'strict': {
        'dir': 'strict',
        'template': 'strict-TEMPLATE.xml',
        'xml': 'strict.xml'
    },
    'logn': {
        'dir': 'logn',
        'template': 'logn-constrained-TEMPLATE.xml',
        'xml': 'logn-constrained.xml'
    },
    'flc': {
        'dir': 'flc',
        'template': 'local-logn-constrained-TEMPLATE.xml',
        'xml': 'local-logn-constrained.xml'
    },
    'local_strict': {
        'dir': 'local-strict',
        'template': 'local-strict-constrained-TEMPLATE.xml',
        'xml': 'local-strict-constrained.xml'
    }
}


def create_xml(sequences, rep_xml, template):
    ff = open(rep_xml, 'w')
    with open(template, 'r') as f:
        for line in f:
            line = line.rstrip('\n').rstrip('\r')
            if 'TEMPLATE_SEQUENCE' in line:
                ff.write(str(sequences))
            else:
                ff.write(line + '\n')
    ff.close()


def create_true_tree(filename):
    #   0.005
    #  -------------E (0.005)
    # |       0.01
    # |       ----- H (0.01)
    # |0.005 |
    #  ------
    #        |
    #         ----- B (0.015)
    #         0.015

    # equine
    E = '(E0_0[&class=0]:2,E1_1[&class=0]:1)'
    # human
    H = '(H0_0[&class=1]:2,H1_1[&class=1]:1)'
    # bird
    B = '(B0_0[&class=2]:2,B1_1[&class=2]:1)'

    for i in xrange(2, 20):
        E = '(' + E + '[&class=0]:1,E' + str(i) + '_' + str(i) + '[&class=0]:1)'
        H = '(' + H + '[&class=1]:1,H' + str(i) + '_' + str(i) + '[&class=1]:1)'
        B = '(' + B + '[&class=2]:1,B' + str(i) + '_' + str(i) + '[&class=2]:1)'

    with open(filename, 'w') as f:
        f.write('#NEXUS\n')
        f.write('begin trees\n')
        f.write('tree TREE1 = [&R] ')
        f.write('(' + E+ '[&class=0]:30,(' + H + '[&class=1]:25,' + B + '[&class=2]:25)[&class=0]:5);\n')
        f.write('end;\n')



def sim(i):
    seed = i
    rep_dir = 'sim{}'.format(i)
    if not os.path.lexists(rep_dir):
        os.mkdir(rep_dir)

    rep_fasta = os.path.join(rep_dir, 'seqs.fa')

    if not os.path.lexists(rep_fasta):
        # Simulate sequences
        rep_tree = os.path.join(rep_dir, 'tree.tree') # contains a record of everything
        cmd = [simultron, '-i', true_tree, '-t', rep_tree, '-o', rep_fasta, '-F', 'fasta', '-m', 'HKY', '-r', '3',
               '-l', '10000', '-d', 'logn', '-M', '0.005,0.01,0.015', '-S', '0.2,0.2,0.2', '-R', str(seed)]
        check_call(cmd)

    sequences = ''
    with open(rep_fasta) as f:
        for record in SeqIO.parse(f, "fasta"):
            sequences += '<sequence id="seq_' + record.id + '" taxon="' + record.id + \
                         '" totalcount="4" value="' + record.seq + '"/>\n'

    for k in data.keys():
        model_dir = os.path.join(rep_dir, data[k]['dir'])
        if not os.path.lexists(model_dir):
            os.makedirs(model_dir)

        # Create BEAST file from template
        create_xml(sequences, os.path.join(rep_dir, data[k]['dir'], data[k]['xml']), data[k]['template'])

        # Run BEAST
        cmd = [os.path.join(beast_dir, 'beast'), '-overwrite', data[k]['xml']]
        check_call(cmd, cwd=model_dir)

        cmd = [os.path.join(beast_dir, 'treeannotator'), '-burnin', str(50), 'data.trees', 'data-annot.tree']
        check_call(cmd, cwd=model_dir)


create_true_tree(true_tree)

inputs = xrange(20, 30)
results = Parallel(n_jobs=3)(delayed(sim)(i) for i in inputs)
