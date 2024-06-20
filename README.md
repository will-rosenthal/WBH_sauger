Scripts and additional figures for "Influence of dams on sauger population structure and hybridization with introduced walleye" 

Below is a guide to the figures included in this repository, including a caption for each.

![ibd_bh_spawn](https://github.com/will-rosenthal/WBH_sauger/assets/39491429/dfbebc37-d5b8-4331-ba9e-5fa1ee00ec68)
Assessment of isolation-by-distance for sauger collected in the Bighorn River during the spawning period. Each point represents a comparison between two sampling locations, with the river distance between them on the x-axis and their divergence (as measured by Fst) on the y-axis. A mantel test to account for spatial autocorrelation of sampling points reveals that this relationship is not significant (p=0.917).

![loc_pca](https://github.com/will-rosenthal/WBH_sauger/assets/39491429/decd4365-ebcb-469b-8ede-352788f62683)
Plots show PCA for a subset of individuals based on their location in the Bighorn River, with each subset labeled above the plot. Each point represents one individual, and each point is colored by the temporal event that individual was sampled during. The percent of total variance explained by each principal component is indicated on the x and y-axis labels. These results show no consistent genetic differentiation between individuals sampled at different times of year within each location.

![pca_all_12_temporalevent_miss0 25](https://github.com/will-rosenthal/WBH_sauger/assets/39491429/ffed1133-e768-4425-8d31-22496990fd68)
There is no consistent genetic differentiation apparent between individuals sampled in different locations during the same time of year. Each plot shows a PCA of a subset of individuals based on the time of year during that individual’s collection, with each subset labeled above the plot. Each point represents one individual, and each point is colored by the location that individual was sampled at. “Non-spawning” signifies an individual that, based on the time of year (typically fall), should not be undergoing a spawning migration. The percent of total variance explained by each principal component is indicated on the x and y-axis labels.

![pca_bighorn_spawn_12_date_miss0 25](https://github.com/will-rosenthal/WBH_sauger/assets/39491429/d0d71a18-af94-445f-9d37-3eb385510cd4)
There is weak clustering of individuals sampled in the same year during the spawning period in the Bighorn River. Each plot shows the same set of points from a PCA of spawning Bighorn River sauger, but each plot colors individuals from a different location. Each point represents one individual, and each point is colored by the date that individual was sampled and shaped according to its sampling location. The percent of total variance explained by each principal component is indicated on the x and y-axis labels.

![pca_wind_12_0 25miss](https://github.com/will-rosenthal/WBH_sauger/assets/39491429/eb860559-e080-497e-9469-e00241e56f75)
There are no consistent genetic differences between individuals sampled at different locations upstream of Boysen Dam. PCA of sauger individuals sampled in the areas upstream of Boysen Dam. Each point represents one individual, with color indicating sample collection location. The percent of total variance explained by each principal component is indicated on the x and y-axis labels.

![pca_23_0 25miss](https://github.com/will-rosenthal/WBH_sauger/assets/39491429/ca910198-f83c-4d4c-ab91-9b61b868c136)
Principal component axes (PCs) 2 and 3 from PCA of all sauger in the dataset are here colored by individual sampling location. Differentiation is evident on PC3, but this does not correspond with sampling location, and the explanation for this structure remains unclear. The percent of total variance explained by each principal component is indicated on the x and y-axis labels.

[DAPC_other_loci-1.pdf](https://github.com/user-attachments/files/15917798/DAPC_other_loci-1.pdf)
The cryptic differentiation evident in PCs 2 and 3 (see above figure) is tied to loci distributed across the genome, not confined to any single region of the genome. Discriminant analysis of principal components results for all sauger. Each individual was binned into either group 1 or 2 depending on its membership to the clustering found on principal components 2 and 3 of a PCA of all sauger samples. Group 1 was all individuals with PC 3 scores less than -0.01, and group 2 was all individuals with PC 3 scores greater than zero. Individuals with PC 3 scores less than zero but greater than -0.01 were excluded. The left window shows the distribution of individuals along the first discriminant function (colored by group membership), and the right window shows the discriminant function loading for each SNP across the genome. The horizontal dashed black line shows the significance threshold obtained from a DAPC randomization procedure. Each SNP with a discriminant function loading higher than the significance threshold is labeled with its location on the scaffold.

![pca_12_lib_miss0 25](https://github.com/will-rosenthal/WBH_sauger/assets/39491429/cd5fbf54-06ba-4305-a493-9ce9ca5d1fcd)
The cryptic differentiation evident in PCs 2 and 3 (Fig \ref{fig:pca3}) does not appear to be tied to sequencing library effects or other potential methodological sources. PCA results for all sauger samples, with each point's color and shape corresponding to the sequencing library that sample was prepared in. Some samples were sequenced in both library 3 and library 4, and are denoted with the label "3+4".

![WBH_diverge_0 3_nomig_nomaf_mod](https://github.com/will-rosenthal/WBH_sauger/assets/39491429/62ad7d85-16fd-484c-9bc4-58158398802a)
Visualization of the residuals from the divergence with no migration demographic models fit in \texttt{dadi}. While this figure only shows the residuals for one round of model fitting, all other rounds that also found the global optimum had similarly structured residuals.
