---
title: "CUSTARD User's Guide"
author: "Kamil Slowikowski"
date: "February 16, 2016"
output: html_document
---



# Introduction

CUSTARD finds pre-defined expression profiles (templates) that are overabundant
in short time course gene expression data.

The algorithm consists of the following steps:

1. Create a small set of expression profile templates.
2. Assign genes to each template.
3. Test each template for enrichment with gene expression profiles.

CUSTARD is an R implementation of the work by Jason Ernst et al. described in:

> Ernst, J., Nau, G. J. & Bar-Joseph, Z. Clustering short time series gene expression data. Bioinformatics 21 Suppl 1, i159–68 (2005). https://www.ncbi.nlm.nih.gov/pubmed/15961453

We have extended the method to meet the specific needs of our analysis, and we
provide those extensions in the CUSTARD R package.

This document is meant to help new users understand the concepts underlying the
method.

# Templates

We want to find overabundant expression profiles in a short time series gene
expression dataset. To do this, we will assume a simple model and create a few
profiles that will serve as templates.

After each time point, we expect that a gene will do one of the following:

- stay the same
- increase in expression
- decrease in expression

In our model, the gene will increase or decrease in discrete units. To better
understand what the landscape of such a model looks like, let's visualize all
possible templates with `n = 5` timepoints and allowing `c = 1, 2, 3` units of
change at each time point:

![plot of chunk all_templates_lines](https://github.com/raychaudhurilab/custard/blob/master/vignettes/figures/all_templates_lines-1.png)

At time point 2, notice that the templates under `c = 1` are allowed to stay
the same as time point 1, increase by one unit, or decrease by one unit. The
gene has 3 choices at each transition between time points, for a total of
$3^4 = 81$ distinct templates. 

We can also visualize the 81 templates for `c = 1` as a heatmap:

![plot of chunk all_templates_heatmap](https://github.com/raychaudhurilab/custard/blob/master/vignettes/figures/all_templates_heatmap-1.png)

Use a greedy algorithm to select `m = 10` distinct profiles from the 81
possibilities when `c = 1`:

![plot of chunk distinct_templates](https://github.com/raychaudhurilab/custard/blob/master/vignettes/figures/distinct_templates-1.png)

Running time of template selection algorithm:

![plot of chunk running_time](https://github.com/raychaudhurilab/custard/blob/master/vignettes/figures/running_time-1.png)

# Enrichment

Test enrichment with 5,000 null genes sampled from the uniform distribution
against 25 templates with `c = 1`:


```r
mat <- matrix(runif(5000 * 5), nrow = 5)
cus <- custard(mat, templates = 25, magnitude = 1)
```

![plot of chunk test_enrichment](https://github.com/raychaudhurilab/custard/blob/master/vignettes/figures/test_enrichment-1.png)

Here is a table of the top 6 templates:


|   |  Expected| Observed|    BinomP|BinomSignif |     PermP|PermSignif |
|:--|---------:|--------:|---------:|:-----------|---------:|:----------|
|8  | 525.81667|      563| 0.0422226|FALSE       | 0.0661157|FALSE      |
|4  |  58.89167|       69| 0.0848910|FALSE       | 0.1157025|FALSE      |
|14 | 137.55000|      152| 0.0995330|FALSE       | 0.1735537|FALSE      |
|34 | 290.00000|      308| 0.1319462|FALSE       | 0.1652893|FALSE      |
|38 | 361.03333|      375| 0.2137760|FALSE       | 0.2727273|FALSE      |
|44 | 148.67500|      157| 0.2294014|FALSE       | 0.2892562|FALSE      |

Let's take a closer look at template 23, the one with the strongest P value:

![plot of chunk enrichment_test_template](https://github.com/raychaudhurilab/custard/blob/master/vignettes/figures/enrichment_test_template-1.png)
