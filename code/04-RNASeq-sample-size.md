RNASeq Sample Size
================
Kathleen Durkin
2023-11-02

The biggest obstacle to calculating desired number of samples for the
cod RNA sequencing is choosing a coefficient of variation ($\sigma$),
which is a parameter we would normally need RNA Seq data to estimate.
The package authors provide estimates of $\sigma$ values for several
taxa, the closest to code being zebrafish. They calculate $\sigma_{90}$
values for each sample (i.e., the value such that 90% of the genes have
a smaller variation). Values for the zebrafish are 0.19 and 0.26, and
human samples range from 0.32 to 0.74 with a median of 0.43 – the
suggested default values in the edgeR user guide of 0.1 and 0.4 for
inbred animal and human studies, in the case where no replicates are
available. (Hart et al. 2013). We’ll calculate necessary sample size for
a range of coefficient of variation values.

We will also use Hart et al.’s recommended depth of 20.

The third important parameter is “effect,” or the minimum change in
expression between groups that we want/will be able to detect. We’ll
test a range of values for detectable expression, to see what our
options are.

For the final two parameters, significance level (alpha) and power,
we’ll use the standard values 0.05 and 0.9.

``` r
rnapower(depth=20, cv=c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4), effect=c(1.25, 1.5, 1.75, 2), alpha=0.05, power=0.9)
```

    ##          1.25       1.5      1.75        2
    ## 0.1  25.32263  7.669561  4.026220 2.624379
    ## 0.15 30.59818  9.267386  4.865016 3.171125
    ## 0.2  37.98394 11.504341  6.039331 3.936568
    ## 0.25 47.47993 14.380426  7.549163 4.920710
    ## 0.3  59.08613 17.895642  9.394514 6.123551
    ## 0.35 72.80256 22.049987 11.575384 7.545089
    ## 0.4  88.62920 26.843463 14.091771 9.185326

This output table shows the sample size needed for a given detectable
expression and coefficient of variation. For example, if we assume cod
will have a coefficient of variation similar to the zebrafish of 0.2 and
we want to detect a 25% change in expression (expression=1.25), we would
need a sample size of 38.

We can also try changing depth to see how that would affect necessary
sample size

``` r
rnapower(depth=100, cv=c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4), effect=c(1.25, 1.5, 1.75, 2), alpha=0.05, power=0.9)
```

    ##           1.25       1.5      1.75        2
    ## 0.1   8.440876  2.556520  1.342073 0.874793
    ## 0.15 13.716424  4.154345  2.180869 1.421539
    ## 0.2  21.102190  6.391301  3.355184 2.186982
    ## 0.25 30.598176  9.267386  4.865016 3.171125
    ## 0.3  42.204381 12.782601  6.710367 4.373965
    ## 0.35 55.920805 16.936947  8.891237 5.795503
    ## 0.4  71.747447 21.730422 11.407625 7.435740

Citations: Hart SN, Therneau TM, Zhang Y, Poland GA, Kocher JP.
Calculating sample size estimates for RNA sequencing data. J Comput
Biol. 2013 Dec;20(12):970-8. doi: 10.1089/cmb.2012.0283. Epub 2013 Aug
20. PMID: 23961961; PMCID: PMC3842884.
