import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import matplotlib.pyplot as plt
from bioinfokit.visuz import cluster

df = pd.read_table("norm_deseq_Thalassiosiralesgenes", sep='\t')
df.head()
df[df.astype(bool).sum(1) > len(df.columns)-5] 
df = df.drop(['Gene'], axis=1)
df.head()
df_standar =  StandardScaler().fit_transform(df)
pd.DataFrame(df_standar, columns=df.columns).head()
pca_out = PCA().fit(df_standar)
pca_out.explained_variance_ratio_
np.cumsum(pca_out.explained_variance_ratio_)
load_pca = pca_out.components_
num_pc = pca_out.n_features_
pc_list = ["PC"+str(i) for i in list(range(1, num_pc+1))]
load_pca_df = pd.DataFrame.from_dict(dict(zip(pc_list, load_pca)))
load_pca_df['variable'] = df.columns.values
load_pca_df = load_pca_df.set_index('variable')
load_pca_df

pca_out.explained_variance_
df1 = load_pca_df[['PC1','PC2']]
df1.reset_index(inplace=True)
zone = pd.read_table("zone.txt", sep='\t')
df2 = pd.merge(df1, zone, on = ['variable'],  how="outer")

plt.figure(figsize=(8,5))
sns.scatterplot(
    x="PC1", y="PC2",
    hue="zone",
    data= df2,
    s=150,
    legend="full",
    palette="colorblind",
    alpha=0.9
)

plt.savefig('PCA_transcriptomes.png')
