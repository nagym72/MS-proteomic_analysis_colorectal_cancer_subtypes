# Statistical analysis of protein counts in 3 subtypes of colorectal tumors 

The protein abundances for three different subtypes of colorectal cancer are given (example3groups.txt) and the question
arose if it would be possible to distinguish the three subtypes based on their relative protein abundances.

A total of 8 patient samples were collected and named a1-a2-a3, b1-b2-b3 and c1-c2 with each letter from a-c corresponding to one distinct subtype.

In order to compare the samples against each other, prior normalization is required.

1) Columns were normalized for the total count of proteins so that after normalization the sum of all values in each colum is equal while the total sum of counts over all samples remains the same.

2) Rows are normalized such that the values have 0 mean and 1 sd for each protein (K-nearest clustering utilized in next step for hierarchical clustering where non normalized large values would dominate the distance metric)

Hierarchical clustering reveals that the relative protein abundances of subtype c cluster together, while for a and b there seems to be no clear case.
a2/a3 and b2/b3 cluster together, yet a1/b1 are mistakenly clustered together.


![alt text](https://github.com/nagym72/MS-proteomic_analysis_colorectal_cancer_subtypes/blob/main/Hierarchical_clustering_normalized_for_all_samples_and_proteins_.jpg)

3) repeating the same procedure but this time only for subtype a and b and those proteins that show significant differences in normalized expression 
 (p-value < 0.01, inverted beta binomial test based on https://doi.org/10.1093/bioinformatics/btp677)


Hierarchical clustering based on only those proteins that were H0 got rejected (p-val 0.01) shows a clear separation between class a and b, hence allowing to fully classify
all 3 subtypes of colorectal cancer based on their protein abundance profiles.

![alt text](https://github.com/nagym72/MS-proteomic_analysis_colorectal_cancer_subtypes/blob/main/Hierarchical_clustering_normalized_for_relevant_p_value_proteins.jpg)
