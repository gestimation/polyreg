
<!-- README.md is generated from README.Rmd. Please edit that file -->

# polyreg

<!-- badges: start -->
<!-- badges: end -->

Direct polytomous regression is a competing risks model that jointly
analyzes multiple competing events and estimates estimate multiplicative
effects of a binary exposure. This model naturally enforces a sum
restriction to cumulative incidence probabilities by reparameterizing
nuisance parameters using polytomous log odds products. Risk ratios,
odds ratios or sub-distribution hazard ratios for each of competing
events are estimated by stratified inverse probability of censoring
weighted estimators under adjustment for covariates in the nuisance and
censoring models.

## Installation

You can install the development version of polyreg from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("gestimation/polyreg")
```

## Usage

Direct polynomial regression can model multiple competing events jointly
and estimates multipilcative effects on cumulative probabilities (or a
risk when there is only one type of event, that is, survival data) from
competing risk data (or survival data) under right-censoring. The models
assumed in polyreg are specified by (1) nuisance.model, (2) cens.model
and (3) exposure, effect measures, and the time point of interest.

1)  nuisance.model specifies the model formula for the nuisance model
    that represents the relationship between outcome and covariates
    other than exposure, in which the outcome should be formatted as
    type Surv or Event that shoyuld include a time variable and a status
    variable. For competing risk data, event codes 1 and 2 are used, and
    for survival data, event code 1 is used, with code 0 indicating a
    censored observation. Another argument outcome.type selects a model
    suitable for the outcome type by selecting COMPETINGRISK and
    SURVIVAL, and must be consistent with the event code of the model
    formula. To obtain valid estimates of the exposure effects,
    covariates in the model format should include confounding factors to
    be adjusted for.

2)  cens.model is a specification for obtaining inverse probability
    weights to adjust for censoring. The same time variable and status
    variable as in (1) should be used for the outcome, but the event
    code to indicate censoring should be 0. If only the intercept is
    specified as the covariate, weights based on the Kaplan-Meier
    estimator are calculated. If a stratified variable is specified
    using strata() as a covariate, dependent censoring is taken into
    account by the stratified Kaplan-Meier estimator.

3)  The effect measure options include risk ratio (RR), odds ratio (OR)
    and sub-distribution hazard ratio (SHR), any of which can be
    selected using effect.measure1 and effect.measure2. The time point
    at which these multiplicative effects are estimated is also
    specified by time.point.

The output of polyreg is a list of coefficient, cov, summary and
summary.full.

## Example

The code below is an example of direct polynomial regression analysis.
Regression coefficients and variance covariance matrix and both exposure
(platelet) and covariates (age and tcell) in regression models for
cumulative incidence probabilities of event 1 and 2 at 24 years are
presented. Dataset bmt is a well-known competing risks dataset that is
included in package mets.

``` r
library(polyreg)
data(bmt)
output <- polyreg(nuisance.model = Event(time, cause)~age+tcell, exposure = 'platelet', cens.model = Event(time,cause==0)~+1, data = bmt, 
effect.measure1='RR', effect.measure2='RR', time.point=24, outcome.type='COMPETINGRISK')
print(output$coefficient)
#> [1] -0.4472068  1.0052315 -1.3167257 -0.4822598 -1.6759005  0.4145390  0.6308143
#> [8]  0.1200117
print(output$cov)
#>               [,1]          [,2]         [,3]          [,4]          [,5]
#> [1,]  0.0200745145  0.0053338619 -0.016367885 -0.0115355705 -0.0006470073
#> [2,]  0.0053338619  0.0179075875 -0.007469460  0.0007063736 -0.0050881169
#> [3,] -0.0163678853 -0.0074694604  0.120004125 -0.0015748967 -0.0057383515
#> [4,] -0.0115355705  0.0007063736 -0.001574897  0.0304237153 -0.0011071326
#> [5,] -0.0006470073 -0.0050881169 -0.005738352 -0.0011071326  0.0066535606
#> [6,]  0.0455363765  0.0513083209 -0.034391995 -0.0092520156 -0.0314900328
#> [7,] -0.0227574688 -0.0147496152  0.037940472  0.0034018555  0.0107875255
#> [8,]  0.0093699814  0.0137060404  0.006532487 -0.0014248493 -0.0166096829
#>              [,6]         [,7]         [,8]
#> [1,]  0.045536376 -0.022757469  0.009369981
#> [2,]  0.051308321 -0.014749615  0.013706040
#> [3,] -0.034391995  0.037940472  0.006532487
#> [4,] -0.009252016  0.003401855 -0.001424849
#> [5,] -0.031490033  0.010787525 -0.016609683
#> [6,]  0.330148842 -0.162984297  0.100965525
#> [7,] -0.162984297  0.162951942 -0.062476954
#> [8,]  0.100965525 -0.062476954  0.052380913
```

The summaries of analysis results in the list of outputs
(e.g. output\$summary.full below) are in accordance with the format of
model summary function. All regression coefficients above are included
in summary.full and model summary may be used to converted to risk
ratios, odds ratios or sub-distribution hazards ratios by selecting
exponentiate=TRUE option. The summaries can be displayed in Viewer with
customized statistics such as p-values or confidence intervals.

``` r
msummary(output$summary.full, statistic = c("conf.int"), exponentiate = TRUE)
```

<table style="width:99%;">
<colgroup>
<col style="width: 27%" />
<col style="width: 35%" />
<col style="width: 35%" />
</colgroup>
<thead>
<tr class="header">
<th></th>
<th>event 1</th>
<th>event 2</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Intercept</td>
<td>0.639</td>
<td>0.187</td>
</tr>
<tr class="even">
<td></td>
<td>[0.484, 0.844]</td>
<td>[0.159, 0.220]</td>
</tr>
<tr class="odd">
<td>age</td>
<td>2.733</td>
<td>1.514</td>
</tr>
<tr class="even">
<td></td>
<td>[2.102, 3.552]</td>
<td>[0.491, 4.668]</td>
</tr>
<tr class="odd">
<td>tcell</td>
<td>0.268</td>
<td>1.879</td>
</tr>
<tr class="even">
<td></td>
<td>[0.136, 0.528]</td>
<td>[0.852, 4.145]</td>
</tr>
<tr class="odd">
<td>platelet</td>
<td>0.617</td>
<td>1.128</td>
</tr>
<tr class="even">
<td></td>
<td>[0.439, 0.869]</td>
<td>[0.720, 1.766]</td>
</tr>
<tr class="odd">
<td>effect.measure</td>
<td>RR of platelet at 24 ( ref = 0 )</td>
<td>RR of platelet at 24 ( ref = 0 )</td>
</tr>
<tr class="even">
<td>n.events</td>
<td>153 events in 408 observations</td>
<td>76 events in 408 observations</td>
</tr>
<tr class="odd">
<td>median.follow.up</td>
<td>7.862 [ 0.03 , 110.625 ]</td>
<td>7.862 [ 0.03 , 110.625 ]</td>
</tr>
<tr class="even">
<td>n.loop.iteration</td>
<td>2</td>
<td><ul>
<li></li>
</ul></td>
</tr>
<tr class="odd">
<td>n.optimization.iteration</td>
<td>0</td>
<td><ul>
<li></li>
</ul></td>
</tr>
<tr class="even">
<td>max.function.value</td>
<td>8.63044214525472e-10</td>
<td><ul>
<li></li>
</ul></td>
</tr>
<tr class="odd">
<td>optimization.message</td>
<td>Function criterion near zero</td>
<td><ul>
<li></li>
</ul></td>
</tr>
</tbody>
</table>

Another example dataset is diabetes.complications. Direct polynomial
regression is applied to this dataset to estimate risk ratios of the
first quartile of fruit intake (fruitq1=1 indicates the first quartile,
fruitq1=0 indicates the third to fourth quartile) on diabetes
complications (diabetic retinopathy coded by epsilon=1 and macrovascular
complications coded by epsilon=2). The code below is survival analysis
(outcome.type=‘SURVIVAL’) to estimate the effects on the risk of
diabetic retinopathy at 8 years of follow-up, treating macrovascular
complications as censored.

``` r
data(diabetes.complications)
diabetes.complications$d <- (diabetes.complications$epsilon>0)
model <- "Event(t,d) ~ age+sex+bmi+hba1c+diabetes_duration+drug_oha+drug_insulin+sbp+ldl+hdl+tg+current_smoker+alcohol_drinker+ltpa"
model <- as.formula(model)
output <- polyreg(nuisance.model = model, exposure = 'fruitq1',
                  cens.model = Event(t,d==0)~strata(strata), data = diabetes.complications,
                  effect.measure1='RR', time.point=8, outcome.type='SURVIVAL')
```

Only the regression coefficient of exposure is included in summary, so
now model summary does not display parameters of the covariates other
than exposure.

``` r
msummary(output$summary, statistic = c("conf.int"), exponentiate = TRUE)
```

|                      | event 1 (no competing risk)              |
|----------------------|------------------------------------------|
| fruitq1 ( ref = 0 )  | 1.366                                    |
|                      | \[1.154, 1.617\]                         |
| effect.measure       | RR of fruitq1 at 8 ( ref = 0 )           |
| n.events             | 358                                      |
| n.events.exposed     | 113 events in 258 exposed observations   |
| n.events.unexposed   | 245 events in 720 unexposed observations |
| median.follow.up     | 7.9973 \[ 0.0493 , 11.0034 \]            |
| n.parameters         | 16                                       |
| optimization.message | Function criterion near zero             |

The code below specifies direct polytomous regression of both of
competing events (outcome.type=‘COMPETINGRISK’). Initial values of
regression parameters are imported as data.initial.values.

``` r
data.initial.values <- c(-20.01, 1.591, -0.2065, 1.962, 10.05, 4.004, 0.9624, 1.582, 2.414, -0.02140,
                         1.072, 1.304, -0.6444, -1.600, -0.01816, 0.4460, -29.25, 6.586, -1.090, 3.515,
                         4.048, -5.012, 0.05348, 0.1609, 6.432, -0.4783, -3.230, -4.135, 1.451, 0.04964,
                         3.822, -0.08824)
model <- "Event(t,epsilon) ~ age+sex+bmi+hba1c+diabetes_duration+drug_oha+drug_insulin+sbp+ldl+hdl+tg+current_smoker+alcohol_drinker+ltpa"
model <- as.formula(model)
output <- polyreg(nuisance.model = model, exposure = 'fruitq1',
                  cens.model = Event(t,epsilon==0)~strata(strata), data = diabetes.complications,
                  effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='COMPETINGRISK', data.initial.values=data.initial.values)
msummary(output$summary, statistic = c("conf.int"), exponentiate = TRUE)
```

<table style="width:99%;">
<colgroup>
<col style="width: 21%" />
<col style="width: 39%" />
<col style="width: 38%" />
</colgroup>
<thead>
<tr class="header">
<th></th>
<th>event1</th>
<th>event2</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>fruitq1 ( ref = 0 )</td>
<td>1.558</td>
<td>0.913</td>
</tr>
<tr class="even">
<td></td>
<td>[1.311, 1.851]</td>
<td>[0.461, 1.806]</td>
</tr>
<tr class="odd">
<td>effect.measure</td>
<td>RR of fruitq1 at 8</td>
<td>RR of fruitq1 at 8</td>
</tr>
<tr class="even">
<td>n.events</td>
<td>279</td>
<td>79</td>
</tr>
<tr class="odd">
<td>n.events.exposed</td>
<td>91 events in 258 exposed observations</td>
<td>22 events in 258 exposed observations</td>
</tr>
<tr class="even">
<td>n.events.unexposed</td>
<td>188 events in 720 unexposed observations</td>
<td>57 events in 720 unexposed observations</td>
</tr>
<tr class="odd">
<td>median.follow.up</td>
<td>7.9973 [ 0.0493 , 11.0034 ]</td>
<td>7.9973 [ 0.0493 , 11.0034 ]</td>
</tr>
<tr class="even">
<td>n.parameters</td>
<td>16</td>
<td>16</td>
</tr>
<tr class="odd">
<td>optimization.message</td>
<td>Function criterion near zero</td>
<td><ul>
<li></li>
</ul></td>
</tr>
</tbody>
</table>
