from re import S
import pandas as pd
pd.plotting.register_matplotlib_converters()
import matplotlib.pyplot as plt # for data visualization
import seaborn as sns # for statistical data visualization

from Bio import SeqIO
from Bio import pairwise2

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

import warnings
warnings.filterwarnings('ignore')

# Reading the sequences files
SARS = SeqIO.read("./src/fasta/sars.fasta", "fasta")
MERS = SeqIO.read("./src/fasta/mers.fasta", "fasta")
COV2 = SeqIO.read("./src/fasta/cov2.fasta", "fasta")
EBOLA = SeqIO.read("./src/fasta/EBOLAV.fasta", "fasta")
HEDGEHOG = SeqIO.read("./src/fasta/hedgehog.fasta", "fasta")
HIV2 = SeqIO.read("./src/fasta/hiv2.fasta", "fasta")
CIVET = SeqIO.read("./src/fasta/Civet-SARS.fasta", "fasta")
ZIKAVIRUS = SeqIO.read("./src/fasta/zikavirus.fasta", "fasta")
MALARIA = SeqIO.read("./src/fasta/plasmodium-malariae.fasta", "fasta")
HCV = SeqIO.read("./src/fasta/HCV.fasta", "fasta")
DEN1 = SeqIO.read("./src/fasta/den1.fasta", "fasta")
DEN2 = SeqIO.read("./src/fasta/den2.fasta", "fasta")
DEN3 = SeqIO.read("./src/fasta/den3.fasta", "fasta")
DEN4 = SeqIO.read("./src/fasta/den4.fasta", "fasta")

# Creating the dataframe
df_data = {
    'ID': [SARS.id, MERS.id, COV2.id, EBOLA.id, HEDGEHOG.id, HIV2.id, CIVET.id, ZIKAVIRUS.id, MALARIA.id, HCV.id, DEN1.id, DEN2.id, DEN3.id, DEN4.id],
    'Name': [SARS.name, MERS.name, COV2.name, EBOLA.name, HEDGEHOG.name, HIV2.name, CIVET.name, ZIKAVIRUS.name, MALARIA.name, HCV.name, DEN1.name, DEN2.name, DEN3.name, DEN4.name],
    'Description': [SARS.description, MERS.description, COV2.description, EBOLA.description, HEDGEHOG.description, HIV2.description, CIVET.description, ZIKAVIRUS.description, MALARIA.description, HCV.description, DEN1.description, DEN2.description, DEN3.description, DEN4.description],
    'Seq': [SARS.seq, MERS.seq, COV2.seq, EBOLA.seq, HEDGEHOG.seq, HIV2.seq, CIVET.seq, ZIKAVIRUS.seq, MALARIA.seq, HCV.seq, DEN1.seq, DEN2.seq, DEN3.seq, DEN4.seq]
}

df_indexs = ['SARS', 'MERS', 'COV2', 'EBOLA', 'HEDGEHOG', 'HIV2', 'CIVET', 'ZIKAVIRUS', 'MALARIA', 'HCV', 'DEN1', 'DEN2', 'DEN3', 'DEN4']
df = pd.DataFrame(df_data, index=df_indexs)

# Selecting the data to create the train dataset and the corresponding dataset to use.
X = df.drop(['Seq'], axis=1)
y = df['Seq']

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state = 0)
print(X_train.shape, X_test.shape)

# ...
cols = X_train.columns
print(cols)
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

X_train = pd.DataFrame(X_train, columns=[cols])
X_test = pd.DataFrame(X_test, columns=[cols])

print(X_train.describe())

#SARS_COV = pairwise2.align.globalxx(SARS.seq, COV2.seq, one_alignment_only=True, score_only=True)
#print('SARS/COV Similarity (%):', SARS_COV / len(SARS.seq) * 100)

#MERS_COV = pairwise2.align.globalxx(MERS.seq, COV2.seq, one_alignment_only=True, score_only=True)
#print('MERS/COV Similarity (%):', MERS_COV / len(MERS.seq) * 100)

#MERS_SARS = pairwise2.align.globalxx(MERS.seq, SARS.seq, one_alignment_only=True, score_only=True)
#print('MERS/SARS Similarity (%):', MERS_SARS / len(SARS.seq) * 100)

#X = ['SARS/COV2', 'MERS/COV2', 'MERS/SARS']
#Y = [SARS_COV/ len(SARS.seq) * 100, MERS_COV/ len(MERS.seq)*100, MERS_SARS/len(SARS.seq)*100]
#plt.title('Sequence identity (%)')
#plt.bar(X,Y,color=(0.2, 0.4, 0.6, 0.6))