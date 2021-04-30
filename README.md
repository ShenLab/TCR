## TCR data analysis and modeling
1. In Data_process, “CDR3 and peptide data preprocess.ipynb” is used to process raw data at the first place. “ERGOdata_process.ipynb”  is used to generate train+validate and test data for ERGO methods. “repetitive data process_ae.ipynb” and “repetitive data process_lstm” are used to remove repetitive data among train, validate and test data used in ERGO methods.

2. ERGO-result contains summaries about data size used in deep learning models and results from those models. ERGO_result_process.ipynb is used to summary results and generate plots based on results.

3. I randomly chose two data sets with different numbers(110,000 and 210,000) from entire data for ERGO methods. “ERGO-master-test10” contains codes and results for the small dataset(110,000), and “ERGO-master-test20” contains codes and results for the dataset(210,000). 

4. In deepTCR, “ERGO_compare_DeepTCR.ipynb” is used to train deepTCR methods, and generate results.

5. In baseline, “data_preprocess.ipynb” is used to prepare data for the baseline method.  “editing disatance.ipynb” is used to acquire results from baseline analyses using edit distance method. 
