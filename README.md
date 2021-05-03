## TCR data analysis and modeling
1. In Data_process, “CDR3 and peptide data preprocess.ipynb” is used to process raw data at the first place. “ERGOdata_process.ipynb”  is used to generate train+validation and test data for ERGO methods. “repetitive data process_ae.ipynb” and “repetitive data process_lstm” are used to remove repetitive data among train, validation and test datasets.

2. ERGO-result contains summaries about data sizes used in ERGO and deepTCR methods and results from them. "ERGO_result_process.ipynb" is used to summary and generate plots based on results.

3. For testing ERGO methods, I randomly chose two data sets with different amounts(110,000 and 210,000) from the entire data. “ERGO-master-test10” contains codes and results for the small dataset(110,000), and “ERGO-master-test20” contains codes and results for the large dataset(210,000). In both files, "ERGO.py" is used to split train and validation data, then create negative and positive examples for them. I use "main" function in "copy_ERGO.py" to train ERGO with train and validation data. "pep_test" function in "new_ERGO.py" is used to test ERGO with the test data, and "pep_test" function in "copy_ERGO.py" is used to test ERGO based on each single peptide in the test data.

4. In deepTCR, “ERGO_compare_DeepTCR.ipynb” is used to generate suitable data format and produce results for deepTCR method.

5. In baseline, “data_preprocess.ipynb” is used to prepare data for the baseline method.  “editing disatance.ipynb” is used to acquire results from baseline analyses using edit distance method. 
