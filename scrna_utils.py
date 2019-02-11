#! /usr/bin/python3.4
import numpy as np

# for writing molecular count matrix with format GID\tSYMBOL\tCTS_CELL1\tCTS_CELL2\t...
def write_matrix(gids,genes,matrix,outfile):
	with open(outfile,'w') as g:
		for i in range(len(gids)):
			vector = matrix[i]
			st = gids[i]+'\t'+genes[i]+'\t'+'\t'.join([str(pt) for pt in vector])+'\n'
			g.write(st)
	return 0


# for loading normalized count matrix (floating values) with format GID\tSYMBOL\tNCTS_CELL1\tNCTS_CELL2\t...
# if srt==1, matrix will be sorted by gid
def load_nmatrix(infile,srt):
	gids = []
	genes = []
	matrix = []
	with open(infile) as f:
		for line in f:
			llist = line.split()
			gids.append(llist[0])
			genes.append(llist[1])
			matrix.append([float(pt) for pt in llist[2::]])
	gids = np.array(gids)
	if srt == 1:
		ind = np.argsort(gids)
		gids = gids[ind]
		genes = np.array(genes)[ind]
		matrix = np.array(matrix)[ind]
	else:
		genes=np.array(genes)
		matrix=np.array(matrix)
	return gids, genes, matrix

# for loading molecular count matrix with format GID\tSYMBOL\tCTS_CELL1\tCTS_CELL2\t...
# if srt==1, matrix will be sorted by gid
def load_matrix(infile,srt):
	gids = []
	genes = []
	matrix = []
	with open(infile) as f:
		for line in f:
			llist = line.split()
			gids.append(llist[0])
			genes.append(llist[1])
			matrix.append([int(pt) for pt in llist[2::]])
	gids = np.array(gids)
	if srt == 1:
		ind = np.argsort(gids)
		gids = gids[ind]
		genes = np.array(genes)[ind]
		matrix = np.array(matrix)[ind]
	else:
		genes=np.array(genes)
		matrix=np.array(matrix)
	return gids, genes, matrix

# for loading molecular count matrix for a list of marker gids with format GID\tSYMBOL\tCTS_CELL1\tCTS_CELL2\t...
# first column in marker_INFILE contains list of marker gids
# if fill==1, counts will be set to zero for any marker absent from the matrix
def load_marker_matrix(matrix_INFILE,marker_INFILE,fill):
	gids = []
	genes = []
	matrix = []
	with open(marker_INFILE) as f:
		markers = set([line.split()[0] for line in f])
	with open(matrix_INFILE) as f:
		for line in f:
			llist = line.split()
			gid = llist[0]
			if gid in markers:
				gids.append(gid)
				genes.append(llist[1])
				try:
					matrix.append([int(pt) for pt in llist[2::]])
				except ValueError:
					matrix.append([float(pt) for pt in llist[2::]])
	if fill == 1:
		for gid in marker_INFILE:
			if gid not in gids:
				gids.append(gid)
				genes.append(gid)
				matrix.append([0 for pt in range(len(matrix[0]))])
	gids = np.array(gids)
	ind = np.argsort(gids)
	gids = gids[ind]
	genes = np.array(genes)[ind]
	matrix = np.array(matrix)[ind]
	return gids, genes, matrix
