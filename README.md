# custard

Find overabundant time series profiles in gene expression data with
permutation testing.

Find pre-defined expression profiles that are overabundant in short time
course gene expression data.

The algorithm consists of the following steps:

1. Create a small set of expression profile templates.
2. Assign genes to each template.
3. Test each template for enrichment with gene expression profiles.

CUSTARD is an R implementation of the work by [Jason Ernst] et al. described
in:

> Ernst, J., Nau, G. J. & Bar-Joseph, Z. Clustering short time series gene
> expression data. Bioinformatics 21 Suppl 1, i159–68 (2005).
> <https://www.ncbi.nlm.nih.gov/pubmed/15961453>

We extended the method to meet the specific needs of our analysis, and we
provide those extensions in the CUSTARD R package.

TODO: Describe those extensions.

[Jason]: http://www.biolchem.ucla.edu/labs/ernst/

