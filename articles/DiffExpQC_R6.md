# Differential Expression Analysis Quality Control.

## Missing Value Analysis

Missing values can reveal potential biases or technical problems in the
data. Figure @ref(fig:missingProtein) summarizes missing protein
abundance estimates per group. Panel A shows how many proteins have
$0 - N$ missing values. Ideally, most proteins are quantified in all
samples within a group. Panel B shows the distribution of mean protein
intensity by the number of missing values. Proteins without missing
values usually have higher average abundance than proteins with one or
more missing values because low-abundance proteins are more likely to
remain undetected. Strongly overlapping distributions can point to other
sources of missingness, such as large sample heterogeneity or technical
problems.

(ref:missingProtein) Missing protein abundance estimates by group. Panel
A: number of proteins with $n$ missing values (`nrNA`). Panel B:
distribution of mean within-group protein intensity by the number of
missing values.

![](DiffExpQC_R6_files/figure-html/missingProtein-1.png)

(ref:missingProtein)

## Variance Analysis of the Data

Panel A in Figure @ref(fig:SDViolin) shows the coefficients of variation
(CVs) for all proteins computed from non-normalized data. Ideally, the
within-group CV should be smaller than the CV across all samples. Panel
B shows the standard deviation distribution for log2-transformed data,
while panel C shows the corresponding distribution after normalization.
Normalization should reduce within-group variance compared with the
overall variance. If normalization increases within-group variance
relative to the overall variance, the selected normalization method may
not be compatible with the data.

(ref:SDViolin) Protein-level variance distributions. Panel A:
coefficients of variation (CVs) within groups and across the full
experiment (`all`) for non-normalized data. Panel B: standard deviations
(SDs) for log2-transformed data. Panel C: SDs after normalization. Black
dots indicate medians.

![](DiffExpQC_R6_files/figure-html/SDViolin-1.png)

(ref:SDViolin)

Table @ref(tab:CVtable) shows the median CV and SD values for all groups
and across all samples (`all`).

| what    |    A |    B | Ctrl |  All |
|:--------|-----:|-----:|-----:|-----:|
| CV      | 3.47 | 3.37 | 3.64 | 9.32 |
| sd_log2 | 0.05 | 0.05 | 0.05 | 0.14 |
| sd      | 0.06 | 0.05 | 0.06 | 0.14 |

Median coefficient of variation (CV) and standard deviation (SD) values.

## Differential Expression Analysis

Most proteins in a dataset are usually not differentially expressed;
therefore, differences between the two groups should be centered close
to zero. Figure @ref(fig:densityOFFoldChanges) shows the distribution of
group differences for all proteins. Ideally, the median of this
distribution (red line) should be close to zero (green line). If the
median and mode of the difference distribution are non-zero, this should
be considered when interpreting the differential expression results.

Panel B in Figure @ref(fig:densityOFFoldChanges) shows the distribution
of p-values for all proteins. If the null hypothesis is true, p-values
should be approximately uniformly distributed. A subset of
differentially expressed proteins produces a higher frequency of small
p-values. A higher frequency of large p-values close to 1 can indicate
that the linear model does not describe the variance of the data well,
for example because of outliers or an unmodeled source of variability.

(ref:densityOFFoldChanges) Differential expression diagnostic
distributions. Panel A: differences between groups for all proteins. The
red dotted line marks the median fold change; the green line marks the
expected median fold change. Panel B: histogram of p-values for all
proteins.

![](DiffExpQC_R6_files/figure-html/densityOFFoldChanges-1.png)

(ref:densityOFFoldChanges)

The MA plot in Figure @ref(fig:MAPlot) helps identify whether large fold
changes are concentrated among high- or low-abundance proteins. Panel A
shows the group difference (y-axis) as a function of average protein
abundance (x-axis). The observed fold change should not depend on
protein abundance. Panel B shows the same group difference against the
rank of average protein abundance.

(ref:MAPlot) MA plot diagnostics. Panel A: group difference versus
average protein abundance. Panel B: group difference versus the rank of
average protein abundance.

![](DiffExpQC_R6_files/figure-html/MAPlot-1.png)

(ref:MAPlot)
