# Tutorial for reproducing the results

1. ERp-newTPM_final.ipynb: our proposed MOPCR model on ERpositive breast cancer
2. HER2E-newTPM_final.ipynb: our proposed MOPCR model on HER2E breast cancer
3. TNBC-newTPM_final.ipynb: our proposed MOPCR model on TNBC breast cancer

For those three codes, first run the first block to do cross-validation and feature selection to get the best model and parameters. （This may take several hours, since it involves the enumeration of all the possible combinations of features and parameters.）

And then it comes to the external validation part. Since all the validation dataset lack some features, please use the code in LASSO folder to do feature imputation.

- LASSO_test_cdkn2a.ipynb: impute methylation of cdkn2a for external datasets
- LASSO_test_e2f.ipynb: impute protein hallmark of e2f targets for external datasets
- LASSO_test_ki67.ipynb: impute ki67 for external datasets
- LASSO_test_kit.ipynb: impute methylation of kit for external datasets
- LASSO_test_map4k1.ipynb: impute methylation of map4k1 for external datasets
- LASSO_test_mtorc1.ipynb: impute protein hallmark of mtorc1 signaling for external datasets
- LASSO_test_uv.ipynb: impute protein hallmark of uv response up for external datasets

After imputation, please come back to the 1-3 ipynb file and run the corresponding block to get the final result.

