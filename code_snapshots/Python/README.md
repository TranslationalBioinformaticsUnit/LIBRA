# Raw code usage:

- LIBRA fine-tuning code ([LIBRA_fine_tune_code](https://github.com/TranslationalBioinformaticsUnit/LIBRA/tree/main/code_snapshots/Python/LIBRA_fine_tune_code))
  - Use paired-omics arrays as input for fine_tune_libra.py
    - Use "fine_tune_libra.py" for LIBRA models training in parallel using desired grid of hyperparameters.
    > **Outputs: Different outputs generated during the training will be stored in the working directory.**

- PPJI metric computation over all fine-tune models:
    - Use **supp_code.R** stored at [R code](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/code_snapshots/R/) employed on LIBRA manuscript.

