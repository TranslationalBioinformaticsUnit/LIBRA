# Raw code usage:

- Preprocessing ([Seurat_code](https://github.com/TranslationalBioinformaticsUnit/LIBRA/tree/main/R/Seurat_code))
  - Removing low-quality features and cells
    - Use "1_pre_analysis_Seurat3.R" (or your own QC pipeline)
    
- Analysis ([Seurat_code](https://github.com/TranslationalBioinformaticsUnit/LIBRA/tree/main/R/Seurat_code))
  - Normalize and visualize with Seurat for independet samples (or other analysis pipeline that generates at least normalized matrix as output). If more than one sample is present for a given omic first integrate them as usual pipelines requires.
    - Use "2_analysis_Seurat3.R" (or your own Normalization pipeline)
    > **Normalization method: Different normalization methods can enhance different LIBRA applications such as - by now - integration or prediction.** As an example Seurat SCT method has shown similar integration performance but a significant decrease on prediction power.

- LIBRA ([Libra_code](https://github.com/TranslationalBioinformaticsUnit/LIBRA/tree/main/R/LIBRA_code))
  - Use analysis_output.RData as input for LIBRA neural network
    - Use "LIBRA.R" for networks training
    - Use "Metrics_LIBRA.R" for additional quality metrix computation
    > **Outputs: Different outputs generated during the training will be stored in the working directory.**

- Use **supp_code.R** for Python fine-tune full list of models PPJI computation employed on LIBRA fine-tune.


