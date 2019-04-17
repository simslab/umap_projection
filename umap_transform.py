#! /usr/bin/python
import sys
import numpy as np
import umap
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] =42
import seaborn as sns

from scipy.stats.stats import spearmanr
from scrna_utils import load_marker_matrix


refrun_NAME = sys.argv[1]      # run name for reference dataset (from which original UMAP embedding was defined)
marker_INFILE = sys.argv[2]    # file containing marker gene list (used to compute original UMAP embedding)
pg_INFILE = sys.argv[3]	       # file containing cluster labels for reference dataset (or any other index)
proj_PREFIX = sys.argv[4]      # prefix for output files containing projection data
k_PARAM = int(sys.argv[5])     # k parameter for UMAP (number of nearest neighbors for knn graph)
used_genes = sys.argv[6]       # output file containing list of marker genes actually used (original list may be filtered)
projrun_NAMES = sys.argv[7::]  # list of run names for projection


print('Loading data...')
rmatrix_INFILE = refrun_NAME+'/'+refrun_NAME+'.matrix.txt'
rgids,rgenes,rmatrix = load_marker_matrix(rmatrix_INFILE,marker_INFILE,1)

pdata = []
for run in projrun_NAMES:
	pmatrix_INFILE = run+'/'+run+'.matrix.txt'
	pdata.append(load_marker_matrix(pmatrix_INFILE,marker_INFILE,1))

print('Filtering data...')
rmatrix_filt = []
pdata_filt = [[] for run in projrun_NAMES]
with open(used_genes,'w') as g:
	for i in range(len(rgenes)):
		if sum(rmatrix[i]) > 0:
			ct=0
			for pdat in pdata:
				if sum(pdat[2][i]) > 0:
					ct+=1
			if ct == len(pdata):
				for j in range(len(projrun_NAMES)):
					pdata_filt[j].append(pdata[j][2][i].tolist())
				rmatrix_filt.append(rmatrix[i])
				gene = rgids[i]
				g.write('%(gene)s\n' % vars())

rmatrix_filt = np.array(rmatrix_filt)
del pdata
del rmatrix

print('Computing model...')
umap_model = umap.UMAP(n_neighbors=k_PARAM,random_state=42,metric='spearman').fit(rmatrix_filt.T)
umap_model_emb = umap_model.embedding_
model_OUTFILE = proj_PREFIX+'.umap_proj_model.txt'
np.savetxt(model_OUTFILE,umap_model_emb,delimiter='\t')

pcolor_set = ['red','green','blue','magenta','brown','cyan','black','orange','grey','darkgreen','yellow','tan','seagreen','fuchsia','gold','olive']
pgs = np.loadtxt(pg_INFILE,dtype='int')
pcolors = np.array([pcolor_set[pg] for pg in pgs])
print('Computing projections...')
for i in range(len(pdata_filt)):
	umap_proj = umap_model.transform(np.float32(np.transpose(np.array(pdata_filt[i]))))
	proj_OUTFILE = proj_PREFIX+'.'+projrun_NAMES[i]+'.umap_proj.txt'
	pdf_OUTFILE = proj_PREFIX+'.'+projrun_NAMES[i]+'.umap_proj.pdf'
	np.savetxt(proj_OUTFILE,umap_proj,delimiter='\t')
	with PdfPages(pdf_OUTFILE) as pdf:
		fig=plt.figure(figsize=(10,10))
		rind = np.argsort(np.random.rand(len(pgs)))  # random index 

		ax0 = fig.add_subplot(2,2,1)
		ax0.scatter(umap_model_emb[rind,0],umap_model_emb[rind,1],c=pcolors[rind],s=1)
		sns.kdeplot(umap_model_emb[rind,0],umap_model_emb[rind,1],cmap='binary_r',shade=False,gridsize=70,n_levels=14)
		ax0.set_aspect('equal')
		ax0.set_axis_off()
		xlm,ylm=ax0.get_xlim(),ax0.get_ylim()
		
		ax1 = fig.add_subplot(2,2,2)
		ax1.scatter(umap_model_emb[rind,0],umap_model_emb[rind,1],c=pcolors[rind],s=1)
		sns.kdeplot(umap_proj[:,0],umap_proj[:,1],cmap='binary_r',shade=False,gridsize=70,n_levels=14)
		ax1.set_xlim(xlm)
		ax1.set_ylim(ylm)
		ax1.set_aspect('equal')
		ax1.set_axis_off()

		ax2 = fig.add_subplot(2,2,3)
		ax2.scatter(umap_model_emb[rind,0],umap_model_emb[rind,1],c='lightgrey',s=1)
		ax2.hexbin(umap_proj[:,0],umap_proj[:,1],cmap='plasma',alpha=0.8,mincnt=1,linewidths=0,edgecolors=None,gridsize=70)
		#plt.colorbar(x,shrink=0.6)
		ax2.set_aspect('equal')
		ax2.set_xlim(xlm)
		ax2.set_ylim(ylm)
		ax2.set_axis_off()
		
		ax3 = fig.add_subplot(2,2,4)
		ax3.scatter(umap_proj[:,0],umap_proj[:,1],c='k',s=1)
		ax3.set_aspect('equal')
		ax3.set_xlim(xlm)
		ax3.set_ylim(ylm)
		ax3.set_axis_off()
		pdf.savefig()
		
