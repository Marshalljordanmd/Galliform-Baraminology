#to put Matthew's kmer analysis program in a loop.
#from 'Matthew_wgks_motif_analysis_rc_k1k2_pcc.py'
#!/usr/bin/python
#I have modified Matthew'a script titled "motif_analysis_rc_k1k2_pcc.py" to eliminate command line inputs
import re, os, sys, getopt
import scipy
from scipy import stats
import os

# if os.path.exists('Galliformes_kmer_matrix.mx'):
# 	os.remove('Galliformes_kmer_matrix.mx')


def revcomp(kmer):
    revc = ""
    revcbp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'M': 'K', 'R': 'Y', 'W': 'S', 'S': 'W', 'K': 'M', 'Y': 'Y'}
    for i in range(len(kmer)):
        b = kmer[i]
        rb = revcbp[b]
        revc = rb + revc
    return revc

def calc_score(motif, a, c, g, t, glen, obs):
    p = 1
    bpp = {'A': a, 'C': c, 'G': g, 'T': t, 'N': 0, 'M': 0, 'R': 0, 'W': 0, 'S': 0, 'K': 0, 'Y': 0, 'B': 0, 'D': 0, 'H': 0, 'V': 0}
    mlen = len(motif)
    bps = list(motif)
    for i in range(mlen):
        bp = bps[i]
        pbp = bpp[bp]
        p = p * pbp
    exp = 2 * glen * p
    score = float((obs - exp))/(obs + exp)
    return score, exp

def intersect(sc1, sc2):
    sc1b = []
    sc2b = []
    for k1 in sc1:
        if k1 in sc2.keys():
            sc1b.append(float(sc1[k1]))
            sc2b.append(float(sc2[k1]))
    return sc1b, sc2b #lists of the scores of kmers common to both genomes

def main():

	samples = ['Anseranas_semipalmata.fasta','Anas_platyrhynchos_whole_genome.fasta','Anas_acuta_whole_genome.fasta','Coturnix_japonica_whole_genome.fasta','Gallus_gallus_whole_genome.fasta','Lagopus_muta.fasta','Phasianus_colchicus.fasta','Tympanuchus_pallidicinctus.fasta','Alectoris_chukar_whole_genome.fasta','Numida_meleagris_whole_genome.fasta','Penelope_pileata_whole_genome.fasta','Cyrtonyx_montezumae.fasta','Callipepla_squamata.fasta','Colinus_virginianus.fasta','Odontophorus_gujanensis.fasta','Alectura_lathami.fasta','Patagioenas_fasciata_whole_genome.fasta','Streptopelia_turtur.fasta','Columba_livia.fasta','Streptopelia_decaocto.fasta']

	#samples = ['rCRS_human.fasta','neander_1.fasta','neander_407.fasta']	
	s_no = int(len(samples))
	kmer_scores = dict()
	results = []
	species = []
			
	kmer = 12
	k_ = kmer
	spec = 'duck/fowl'
	
	for s in range(s_no):	
		print('computing kmer scores for ',samples[s])
		inputfile = samples[s]
		scores_kmer = dict()				
		pos = 0
		n = dict()
		bplist = list()
		acgt = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'M': 0, 'R': 0, 'W': 0, 'S': 0, 'K': 0, 'Y': 0, 'B': 0, 'D': 0, 'H': 0, 'V': 0}
		genome_len = 0
		nchr = 0
		
		fi = open(inputfile,'r')
		for line in fi:
			lin = line.rstrip('\n')
			line = str(lin)
			x = re.search(">",line)
			if x :
				pos = 0; kmer = ""; bplist = list(); nchr = nchr + 1
			else:
				w = list(line)
				for b in w:
					acgt[b.upper()] = acgt[b.upper()] + 1
					pos = pos+1
					genome_len = genome_len + 1
	
					bplist.append(b)
					if pos >= k_ + 1:
						bplist.pop(0)
	
					if pos >= k_:
						kmer = ''.join(bplist).upper()
						rc_kmer = revcomp(kmer)
						z = re.search("[NMRWSYKBDHV]",kmer)
						if not z:
							if kmer in n:
								n[kmer] = n[kmer] + 1
							else:
								n[kmer] = 1
							if rc_kmer in n:
								n[rc_kmer] = n[rc_kmer] + 1
							else:
								n[rc_kmer] = 1                            
		fi.close()
			
		tacgt = acgt['A'] + acgt['C'] + acgt['G'] + acgt['T']
		ap = float(acgt['A'])/tacgt
		cp = float(acgt['C'])/tacgt
		gp = float(acgt['G'])/tacgt
		tp = float(acgt['T'])/tacgt
		
		for km in sorted(n):
			score, exp = calc_score(km, ap, cp, gp, tp, genome_len, n[km])
			scores_kmer[km] = score

		file_list = samples[s].split('_')
		name = file_list[0]+ '_' +file_list[1]
		species.append(name)
		results.append(scores_kmer)

	
	#now to compute PCC for the pairs of sample kmer scores
	
	outputfile = ''
	
	fo = open(outputfile,'a')	
	for sp in range(s_no):
		fo.write(species[sp])
		if sp < s_no-1:
			fo.write('\t')
	
	
	for i in range(s_no):
		fo.write('\n')
		for ii in range(s_no):
			scores1 = results[i] #because each element of the results list is a dictionary, scores1 is a dictionary of keys=kmer, values=score
			scores2 = results[ii]			
			(scores1b, scores2b) = intersect(scores1, scores2) #returns 2 lists containing kmer scores of kmers common to both genomes
			(r, p) = scipy.stats.pearsonr(scores1b, scores2b) #r is Pearson coefficient, p = probability
			r = round(r,4)
			fo.write(str(r))
			if ii < s_no - 1:
				fo.write('\t')
			
	fo.close()

if __name__ == "__main__":
    main()





