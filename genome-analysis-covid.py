import pandas as pd
pd.plotting.register_matplotlib_converters()
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio import pairwise2

SARS = SeqIO.read("./src/fasta/sars.fasta", "fasta")
MERS = SeqIO.read("./src/fasta/mers.fasta", "fasta")
COV2 = SeqIO.read("./src/fasta/cov2.fasta", "fasta")

SARS_COV = pairwise2.align.globalxx(SARS.seq, COV2.seq, one_alignment_only=True, score_only=True)
print('SARS/COV Similarity (%):', SARS_COV / len(SARS.seq) * 100)

MERS_COV = pairwise2.align.globalxx(MERS.seq, COV2.seq, one_alignment_only=True, score_only=True)
print('MERS/COV Similarity (%):', MERS_COV / len(MERS.seq) * 100)

MERS_SARS = pairwise2.align.globalxx(MERS.seq, SARS.seq, one_alignment_only=True, score_only=True)
print('MERS/SARS Similarity (%):', MERS_SARS / len(SARS.seq) * 100)

X = ['SARS/COV2', 'MERS/COV2', 'MERS/SARS']
Y = [SARS_COV/ len(SARS.seq) * 100, MERS_COV/ len(MERS.seq)*100, MERS_SARS/len(SARS.seq)*100]
plt.title('Sequence identity (%)')
plt.bar(X,Y,color=(0.2, 0.4, 0.6, 0.6))