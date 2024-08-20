# Bayesian Composite Likelihood method with the open-faced sandwich adjustment (OFS) for Multivairate network meta-analysis

This repo is to implement Bayesian Composite Likelihood method with open-faced adjustment (BCL-OFS) in our paper ["Exploiting multivariate network meta-analysis: A calibrated Bayesian composite likelihood inference" by Yifei Wang, Lifeng Lin and Yulun Liu](https://www.medrxiv.org/content/10.1101/2024.06.25.24309477v1) 

The proposed method eliminated the need to specify a full likelihood function while allowing for the unavailability of within-study correlations among treatments and outcomes. Additionally, we developed a hybrid Gibbs sampler algorithm along with the Open-Faced Sandwich post-sampling adjustment to enable robust posterior inference.

The repo implemnted the proposed method on a real-world dataset comparing treatment effect of root coverage procedures (Buti 2013). To get the code works, please download all the files under the same folder and set the working directory to the source file location. Then,

1. Running **netowrk_plot.R** to plot the networks
3. Running **mcmc.R** to generate posterior samples and perform OFS adjustment
4. Running **concordance_plot.R** to generate concordance plots
5. Running **sucra_plot.R** for generating SUCRA ranking plots

Other files are used for

1. Main functions for **mcmc.R**: **gibbs_within_mh_equal_tau.R**: 
2. Data files: _CALgain.csv_; _CRC.csv_; _KTgain.csv_; _RecRed.csv_; _trt_mapping.csv_; _clean_data.csv_
3. Results from the previous study (Buti 2013): _paper_estimate.csv_

For the result interpreation and the detailed information, you can find in [our paper](https://www.medrxiv.org/content/10.1101/2024.06.25.24309477v1).

Reference

Buti J, Baccini M, Nieri M, La Marca M, Pini-Prato GP. Bayesian network meta-analysis of root coverage procedures: ranking efficacy and identification of best treatment. J Clin Periodontol. 2013 Apr;40(4):372-86. doi: 10.1111/jcpe.12028. Epub 2013 Jan 24. PMID: 23346965.