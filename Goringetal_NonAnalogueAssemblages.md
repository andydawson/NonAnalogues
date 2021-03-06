No-Analogues in North American Pollen Space
========================================================

Abstract
------------------------
Local vs. regional turnover.  Transient novelty? Temporally adjacent novelty (TAN)?

Introduction
------------------------
Our investigation of non-analogue vegetation, both in the past and in the future is tied to differences from modern conditions (Williams and Jackson, 2007).  A major concern with no-analogue climate in the future is the development of forest communities for which we have no modern analogues as tree species respond individualistically to changing climatic conditions. No-analogues provide a challenge for ecologists, land managers and conservationists, tied to our lack of understanding, both in how the process of change will unfold, and in how ecosystem services are provided by these vegetation communities for which we have no analogue.

Estimates of dissimilarity and turnover in paleoecological records have been used to give us an indication of the novelty of ecosystems with respect to modern vegetation, or to provide information about how much change occurs through time within a record (Williams et al., 2001), often relating these changes to climate.  The challenge in this case is that, as with modern ecosystems, changes in pollen composition at one location may represent points along a successional trajectory, where each time point represents significant compositional turnover, but, relative to the surrounding, heterogeneous landscape, very little change in forest composition occurs.  Thus turnover as a single metric can provide results that are not uniquely tied to overall, landscape scale, changes in vegetation.

In a similar way, dissimilarity of past pollen assemblages with respect to modern pollen (d<sub>m</sub>) cannot indicate changes from one community type to another through time, only dissimilarity from modern ecosystems.  While broad scale spatial changes may be visible, a site that transitions from prairie to forest can retain a similarly low d<sub>m</sub> value, provided near analogues exist in the modern pollen data.  A significant challenge is the fact that in many cases modern pollen reflects an anthropogenic alteration of the landscape relative to fossil assemblages.  The presence of taxa such as Artemisa point to significant opening of woodland canopies for example [REF].

A fuller understanding of vegetation change, and of successional processes may be derived from the use of multiple metrics for vegetation change using pollen data. Three dissimilarity measures can be used to help us understand changes through time:
1.  Dissimilarity from modern.
2.  Turnover with an individual core.
3.  Turnover with respect to the surrounding landscape.

Dissimilarity from modern has been widely used in the paleoecological literature, it represent a measure of how different past ecosystems were to modern.  High dissimilarities are most common in North America at ~14kyr, when it is presumed that a combination of climatic and biotic factors gave rise to new vegetation communities, communities that have no modern analogue.

Turnover within a core just tells you how different one time period is to the next.

Landscape turnover has not been previously used.  It is a measure of how different an assemblage at one time period is to all pollen assemblages within a previous time period.  In this sense, it is an estimate of how novel an ecosystem is with respect to all ecosystems that previously occupied the landscape.  In a sense it is d<sub>m</sub> turned on its head.

To assist in this endevour, we examine pollen records from the Neotoma Database, searching through the past, from the late-glacial to the modern to examine in which cases pollen assemblages appear to be non-analogue from the previous time period.  The conceptual framework suggests that we were looking at landscape level turnover, not site turnover.  While site-level turnover is a commonly used metric in paleoecological analysis, it is hampered by its sensitivity to stochastic effects.  Large scale changes in vegetation at a single site might not reflect the kinds of broad scale community change we expect to see under conditions of future climate change.  As such, comparing dissimilarity at a site to changes at a number of sites across the continent provides a better metric of landscape-scale dissimilarity.

![Alt text](/HandMadeFigures/ConceptualPlot.png)
**Figure 1**. *Turnover estimates for time t_i are measured as the minimum within site dissimilarity from the prior timestep (d_si) or the minimum dissimilarity from the landscape of pollen assemblages at the earlier timestep (d_li).  In all cases the earlier timestep is defined as a time preiod between 250 and 750 years prior to the sample of interest.  These two values then form a continuum of possibilities.  High ds and low dl indicate disturbance events within a stable, heterogeneous landscape, high dl & high ds indicate broad-scale, rapid landscape turnover, high dl, low ds is relatively infrequent, but may indicate local refugia within a shifting landscape*

We expect that ecosystems will show patterns of change that can be explained through the multivariate use of these dissimilarity measures.  For example, rapid species migration during the late-Pleistocene/eraly-Holocene transition should be apparent as high local turnover partnered with relatively low landscape change if entire communities are migrating northward.  Strongly individualistic species responses would be visible as paired high site and landscape turnovers.  Modern successional changes should be visible as sets of high self, low landscape changes.  Low landscape and low self changes should be static ecosystems.

Methods
------------------------
We compile records of pollen from depositional records in the Neotoma Database.  To assess whether pollen assemblages are 'non-analogue' we estimate squared-chord distance from pollen samples to a reference set that includes (1) all samples in the Neotoma Database that are between 250 and 750 calibrated years older than the sample, and (2) includes samples from the reference site.




Pollen data from Neotoma was accessed on Nov 09, 2013 using the `neotoma` pacakge for R (Goring, 2013; `http://www.github.com/ropensci/neotoma`).  The dataset includes 560 sites from across eastern North America (east of 100^o W; Figure 1a), with 21397 samples younger than 21kyr cal. BP (Figure 1b).

![plot of chunk Figure1Plots](figure/Figure1Plots.png) 

**Figure 2**. *Sample plot locations and bin sizes for each age class*.

Because sample size may affect our ability to calculate the 95% CI we also use the squared-chord dissimilarity estimate reported in Gill et al. (2009) of XXX as a secondary check.  This allows us to detect no-analogues using multiple methods.

To determine dissimilarity ofver time we estimate dissimilarity from the data using a bootstrap approach for which a sample is compared against a 'landscape' of sites that are between 250 and 750 years older than the sample in question.  For any sample we first test whether a sample fom the site exists between 250 - 750 prior to the sample of interest.  This will prevent anomalously high analogue distances for sites that have never previously been sampled, particularly when they represent new ecoregions.  For each acceptable site we sample one assemblage from each site with samples in the previous 250 - 750 years.  This produces a single sample from each site in the previous time window from which we estimate the minimum suqared chord dissimilarity.  Since some sites have multiple samples in each 500 year time window we re-sample (with replacement) 100 times for each focal site, producing a sample of 100 minimum (squared-chord) dissimilarity values for each pollen assemblage at each site, for which there is a prior sample.





Results
-------------------------
Pollen samples:
Of the 560 pollen sites obtained from Neotoma for this analysis, 560 sites had assemblages that met our criteria.  For these sites there were 18850 unique assemblages spanning the last 21kyr, approximately 90% of the total assemblages for the sites that met our criteria.  Samples excluded from analysis occur throughout the record

Self and Landscape Dissimilarity:

```
## Error: object of type 'closure' is not subsettable
```

```
## Error: ggplot2 doesn't know how to deal with data of class function
```

**Figure X**. *Turnover and dissimilarity for individual pollen records from the Neotoma database, shows a relationship between the two variables, but also considerable noise.*

Turnover within the core is correlated to landscape dissimilarity (r<sub>s</sub> = 

```

Error in p(cor.est$estimate) : object 'cor.est' not found

```

, p < 0.001), but shows a great deal of variability around the 1:1 line (Figure X).  
The analysis produces a somewhat surprising result.  While dissimilarity is high at the beginning of the Holocene, the most rapid rise in turnover occurs in the modern era, even though the density of sites is higher at this time.  High turnovers are seen at 10kyr, between 6 and 7 kyr and then again in the modern period.  While the no-analogue period of the late-glacial has high dissimilarity in relation to modern time, the actual turnover is not significantly higher than during the Holocene transition.


![plot of chunk dissVsAge](figure/dissVsAge.png) 

**Figure 3**. *Turnover through time in the Neotoma database.*


```r

# This is some basic hypothesis testing.  We're interested in knowing how
# much these change.
model0 <- gam(log(self, 10) ~ 1, family = gaussian, data = rf[rf$age < 15000 & 
    rf$age > -50, ])
model1 <- gam(log(self, 10) ~ s(age), family = gaussian, data = rf[rf$age < 
    15000 & rf$age > -50, ])
model1.500 <- gam(log(self, 10) ~ s(age, k = 40), family = gaussian, data = rf[rf$age < 
    15000 & rf$age > -50, ])

model2.500 <- gam(log(self, 10) ~ s(age, k = 40, by = factor(domain)), family = gaussian, 
    data = rf[rf$age < 15000 & rf$age > -50, ])

model3.500 <- gam(log(self, 10) ~ s(age, k = 40, by = factor(division)), family = gaussian, 
    data = rf[rf$age < 15000 & rf$age > -50, ])
```


Discussion
---------------------------
Turnover or analogues?  When we are looking at non-analogues it turns out that it's really the modern.
