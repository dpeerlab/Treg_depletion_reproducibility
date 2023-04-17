import pandas as pd
import scipy
import schpf

outbase = '/data/peer/sharmar1/factors_tcells/experiment2/output/endo/iter_20/'

inbase = '/data/peer/sharmar1/factors_tcells/experiment2/data/endo/'

endo_coo = scipy.sparse.load_npz(inbase + 'endo_coo_matrix.npz')
gene_names = pd.read_csv(inbase + 'endo_gene_names.csv', index_col = 0).index
cell_names = pd.read_csv(inbase + 'endo_cell_names.csv', index_col = 0).index

for nf in range(14, 28, 2):
    print('Computing ' + str(nf) + ' factors', flush = True)
    for niter in range(10):
        print('Iteration ' + str(niter), flush = True)
        endo_hpf = schpf.scHPF(nfactors=nf, verbose=True)
        endo_hpf.fit(endo_coo)

        pd.DataFrame(endo_hpf.gene_score(), index = gene_names).to_csv(outbase + 'gs_' + str(nf) + 'f_iter_' + str(niter) + '_endo_v1.csv')
        pd.DataFrame(endo_hpf.cell_score(), index = cell_names).to_csv(outbase + 'cs_' + str(nf) + 'f_iter_' + str(niter) + '_endo_v1.csv')
        
        del endo_hpf
