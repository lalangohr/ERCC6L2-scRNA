"""
Part 4 of scRNA-seq data analysis. 

Data integration using scVI and cell type annotation using scANVI.
"""

import scvi
import os
import torch
import matplotlib.pyplot as plt

print(f'scVI version: {scvi.__version__}')
print(f'torch version: {torch.__version__}')


# Folder in which result is written
dir_results = './results'


def main():
        
    # Set up GPU
    use_gpu = True
    device = 'cpu'
    if torch.cuda.is_available():
            print(f"CUDA version: {torch.version.cuda}")
            print(torch.zeros(1).cuda())
            cuda_id = torch.cuda.current_device()
            print(f"ID of current CUDA device: {cuda_id}")
            device = 'cuda:' + str(cuda_id)
            print(f"Name of current CUDA device: {torch.cuda.get_device_name(cuda_id)}")
    else:
            print('CUDA not available')
            use_gpu = False


    # Parse adata
    adata_file = os.path.join(dir_results, 'adata_filtered_seeds.h5ad')
    adata = scvi.data.read_h5ad(adata_file)


    # Specify parameters
    n_hidden=128
    n_layers=2
    n_latent=30
    max_epochs_scVI = 500
    max_epochs_scANVI = 200

    plan_kwargs = {
        "reduce_lr_on_plateau": True,
        "lr_patience": 8,
        "lr_factor": 0.1,
    }

    # Single-cell Variational Inference (scVI), see Lopez et al. 2018

    # Set up the AnnData object for scVI
    # None of the data in adata are modified. Only adds fields to adata.
    scvi.model.SCVI.setup_anndata(adata,
                                labels_key="seed_labels",
                                batch_key="sample",
                                categorical_covariate_keys=["protocol", "year"]) 

    scvi_model = scvi.model.SCVI(adata,
                                n_hidden=n_hidden,
                                n_latent=n_latent,
                                n_layers=n_layers,
                                gene_likelihood="nb", # zinb, nb, poisson, default: zinb
                                encode_covariates=True,
                                deeply_inject_covariates=False,
                                use_layer_norm='both',
                                use_batch_norm='none')

    # Move model to device
    scvi_model.to_device(device)

    # Train the scVI model
    scvi_model.train(max_epochs=max_epochs_scVI,
                    plan_kwargs=plan_kwargs,
                    early_stopping=True,
                    early_stopping_monitor='elbo_validation',
                    early_stopping_patience=10,
                    early_stopping_min_delta=0.0,
                    check_val_every_n_epoch=1,
                    use_gpu=use_gpu)

    fig, axs = plt.subplots(2)
    elbo = scvi_model.history["elbo_validation"]
    elbo.to_csv(os.path.join(dir_results, "elbo_val.csv"))
    axs[0].plot(elbo, label="elbo")
    axs[0].legend()


    # Train scANVI and transfer the labels.

    # Initialize scanVI model with weights from pretrained scVI model.
    scanvi_model = scvi.model.SCANVI.from_scvi_model(scvi_model, 'Unknown')

    # Train the scANVI model
    scanvi_model.train(max_epochs=max_epochs_scANVI,
                    early_stopping=True,
                    early_stopping_monitor="classification_loss_train",
                    train_size=1.0,
                    early_stopping_patience=10,
                    early_stopping_min_delta=0.001,
                    check_val_every_n_epoch=1,
                    use_gpu=use_gpu)

    clt = scanvi_model.history["classification_loss_train"]
    clt.to_csv(os.path.join(dir_results, "class_loss_train.csv"))
    axs[1].plot(clt, label="class. loss train")
    axs[1].legend()
    fig.savefig(os.path.join(dir_results, "training.png"))

    # Predict the missing cell types
    adata.obs["C_scANVI"] = scanvi_model.predict(adata)

    # Get inferred cell type probabilities for each cell
    df_probs = scanvi_model.predict(adata, soft=True)
    n_celltypes=45
    for i in range(1,n_celltypes+1):
            adata.obs["probs_scANVI_" + str(i)] = df_probs[str(i)]

    # Get the latent space
    adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata)

    # Change cell types back into long format
    cts = ['HSCs & MPPs', 'Erythro-myeloid progenitors (EMP)', 'Early erythroid progenitor', 'Late erythroid progenitor',
    'Aberrant erythroid cells', 'Megakaryocyte progenitor (MkP)', 'Eosinophil/Basophil progenitor (EoBaso)', 
    'Lympho-myeloid prog', 'Early promyelocytes', 'Late promyelocytes', 'Myelocytes', 'Classical Monocytes',
    'Non-classical monocytes', 'Plasmacytoid dendritic cell progenitors', 'Plasmacytoid dendritic cells (pDC)',
    'Conventional dendritic cells 1 (cDC1)', 'Conventional dendritic cells 2 (cDC2)',
    'Pre-pro B cells', 'Cycling pro-B and pre-B cells', 'Non-cycling pro-B and pre-B cells', 
    'Small pre-B cells (light chain re-arrangement)', 'Immature B cells', 'Mature naïve B cells',
    'Non-switched memory B cells', 'Class-switched memory B cells', 'CD11c+ memory B cells',
    'Plasma cells',
    'CD4+ naïve T cells', 'CD4+ memory T cells', 'CD4+ cytotoxic T cells', 'CD4+ CD69+ T cells',
    'CD8+ naïve T cells', 'CD8+ tissue-resident memory T cell', 'CD8+ central memory T cells ', 
    'CD8+ effector memory T cells', 'gd-T cells',
    'Natural killer T cells', 'CD56dim CD16+ natural killer cells', 'CD56bright CD16- natural killer cells',
    'Mesenchymal cells 1', 'Mesenchymal cells 2',
    'Natural killer cell progenitor',
    'Immature-like blasts', 'Dendritic-cell like blast', 'Monocyte-like blast']
    dict_Triana = {ct: str(idx + 1) for idx, ct in enumerate(cts)}
    dict_inv = {v: k for k, v in dict_Triana.items()}
    
    # Add cell type to adata
    adata.obs["celltype"] = adata.obs["C_scANVI"].map(dict_inv)


    # Print overview of the model
    print('Overview of scVI model:')
    print(scvi_model)
    print('Overview of scANVI model:')
    print(scanvi_model)


    # Save models
    scvi_model.save(os.path.join(dir_results, 'scvi_model/'), overwrite=True)
    scanvi_model.save(os.path.join(dir_results, 'scanvi_model/'), overwrite=True)

    # Save result
    adata.write(os.path.join(dir_results, 'adata_filtered_seeds_scANVI.h5ad'))

    print('Finished')


if __name__ == "__main__":
    main()
