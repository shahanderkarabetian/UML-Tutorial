# Unsupervised machine learning tutorial materials
Tutorial materials associated with the lecture and demonstrations video on unsupervised machine learning in species delimitation.

Video link: https://www.youtube.com/watch?v=MDVT5xffCNo

The paper: https://doi.org/10.1016/j.ympev.2019.106562


Four files are included:

- PCA-DAPC-RF-tSNE_NMNH_workshop.R : This is the R script for running PCA-DAPC, Random Forest, and t-SNE.
- Metano_UCE_SNPs_forRscript.str : This is the input dataset associated with the R script. SNP data.

The R scripts were tested and used with these R package versions:
adegenet 2.1.1 
randomForest 4.6-14
tsne 0.1-3
PCDimension 1.1.9
mclust 5.4.1
cluster 2.0.7-1
factoextra 1.0.5

- VAE_NMNH_workshop.ipynb : The script for running VAE. Alternate collab link here:https://colab.research.google.com/drive/1ynXdWKzhgcdCPkDFhyttqs_SOdBfaFa8?usp=sharing
- Metano_UCE_SNPs_70percent_onehot.txt : This is the input dataset for use with the above VAE script. One-hot encoded SNP data.

This script will be run through https://colab.research.google.com/. Once you navigate to the page, go to the “Upload” tab and open (or drag/drop) the file, and it should then open a new window with the script. If you do not have or want a google account, I have also run it through Kaggle with a minor modification to the script. I can send you the modified script.
