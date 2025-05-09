# CIMA-CLM model

The CIMA-CLM model is cell type-specific and begins with two parallel encoder branches. The first branch encodes input chromatin sequences into DNA embeddings using the pretrained HyenaDNA model, while the second branch processes single-cell gene expression data to generate RNA embeddings using the pretrained scGPT model (see Methods). The weights of these pretrained models are frozen with stop-gradient to speed up model training. To improve the model's ability to generate more informative embeddings, two additional transformer-based encoders further process the DNA and RNA embeddings, respectively. CIMA-CLM incorporates a fusion decoder that integrates the DNA and RNA embeddings using multi-head cross-attention layers. This allows the model to capture relationships between chromatin sequences and gene expression data. Ultimately, the model predicts chromatin accessibility, measured by chromatin peaks, based on the fused features from both data modalities. 


