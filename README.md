
<!-- README.md is generated from README.Rmd. Please edit that file -->

# polyreg

<!-- badges: start -->
<!-- badges: end -->

The polyreg package implements direct polynomial regression, a model
that estimates multiplicative effects (e.g., risk ratio, odds ratio, or
sub-distribution hazard ratio) on risks, cumulative incidence
probabilities at a specific time point, or common effects over time. A
key feature of this model is its ability to analyze multiple competing
events simultaneously while ensuring that the probabilities sum to one.
This is achieved by reparameterizing nuisance parameters using
polytomous log odds products. Additionally, the package supports direct
binomial regression for survival outcomes and the Richardson model for
binomial outcomes, both of which use log odds products.

The models in polyreg are specified by three main components:

(1)Nuisance model: Describes the relationship between outcomes and
covariates (excluding exposure).

(2)Effect measures and time points: Defines the exposure effect to be
estimated and the time point of interest.

(3)Censoring adjustment: Specifies strata for inverse probability
weighting to adjust for dependent censoring.

## Model specification

### 1. Nuisance Model

The nuisance.model argument specifies the formula linking the outcome to
covariates. Its format depends on the outcome type:

(1)Competing risks or survival outcome: Use Surv() or Event() with time
and status variables.

(2)Binomial outcome: Use standard R formula notation.

Default event codes:

(1)Competing risks outcome: 1 and 2 for event types, 0 for censored
observations.

(2)Survival outcome: 1 for events, 0 for censored observations.

(3)Binomial outcome: 0 and 1.

Event codes can be customized using code.event1, code.event2, and
code.censoring. The outcome.type argument must be set to:

(1)Effects on cumulative incidence probabilities at a specific time:
‘COMPETINGRISK’

(2)Effects on a risk at a specific time: ‘SURVIVAL’

(3)Effects on a risk of a binomial outcome: ‘BINOMIAL’

(4)Common effects on cumulative incidence probabilities over time:
‘PROPORTIONAL’

Covariates included in nuisance.model should adjust for confounding
factors to obtain unbiased exposure effect estimates.

### 2. Effect measures and time points

Three effect measures available:

(1)Risk Ratio (RR)

(2)Odds Ratio (OR)

(3)Sub-distribution Hazard Ratio (SHR)

Set the desired measure using effect.measure1 and, for competing risks
analysis, effect.measure2. The time.point argument specifies the
follow-up time at which effects are estimated.

### 3. Censoring adjustment

Inverse probability weights adjust for dependent censoring. Use the
strata argument to specify stratification variables. If no strata are
specified, Kaplan-Meier weights are used.

## Output

The main components of the output list include:

(1)coefficient: Regression coefficients

(2)cov: Variance-covariance matrix

(3)diagnosis.statistics: Inverse probability weights, influence
functions, and predicted values

(4)summary: Summary of estimated exposure effects

Use the summary output with msummary() to display formatted results. The
regression coefficients and their variance-covariance matrix are
provided as coefficient and cov, respectively, with the first element
corresponding to the intercept term, subsequent elements to the
covariates in nuisance.model, and the last element to exposure. Finally,
diagnosis.statistics is a dataset containing inverse probability
weights, influence functions, and predicted values of the potential
outcomes of individual observations.

## Example 1. Unadjusted competing risks analysis

For the initial illustration, unadjusted analysis focusing on cumulative
incidence probabilities of event 1 and 2 at 8 years is demonstrated.
Regression coefficients and variance covariance matrix of both exposure
(fruitq1) and covariates (intercept in this case) in the fitted direct
polytomous regression are presented.

``` r
library(polyreg)
data(diabetes.complications)
output <- polyreg(nuisance.model = Event(t,epsilon) ~ 1, exposure = 'fruitq1', data = diabetes.complications,
          effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='C', report.nuisance.parameter = TRUE)
#> 
#> Attaching package: 'boot'
#> The following object is masked from 'package:survival':
#> 
#>     aml
print(output$coefficient)
#> [1] -1.38313105  0.30043925 -3.99147261  0.07582589
print(output$cov)
#>              [,1]          [,2]          [,3]          [,4]
#> [1,]  0.007240584 -4.543030e-03 -1.051296e-03  0.0050890586
#> [2,] -0.004543030  9.669548e-03 -1.309544e-05 -0.0005223235
#> [3,] -0.001051296 -1.309544e-05  2.436955e-02 -0.0838268451
#> [4,]  0.005089059 -5.223235e-04 -8.382685e-02  0.2889784118
```

The summaries of analysis results in the list of outputs
(e.g. output\$summary below) are in accordance with the format of model
summary function. All regression coefficients above are included in
summary by setting report.nuisance.parameter = TRUE. Model summary may
be used to converted to risk ratios, odds ratios or sub-distribution
hazards ratios using exponentiate=TRUE option. The summaries can be
displayed in Viewer with customized statistics such as p-values or
confidence intervals.

``` r
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
<td>Intercept</td>
<td>0.251</td>
<td>0.018</td>
</tr>
<tr class="even">
<td></td>
<td>[0.212, 0.296]</td>
<td>[0.014, 0.025]</td>
</tr>
<tr class="odd">
<td>fruitq1 ( ref = 0 )</td>
<td>1.350</td>
<td>1.079</td>
</tr>
<tr class="even">
<td></td>
<td>[1.114, 1.637]</td>
<td>[0.376, 3.094]</td>
</tr>
<tr class="odd">
<td>effect.measure</td>
<td>RR of fruitq1 at 8 ( ref = 0 )</td>
<td>RR of fruitq1 at 8 ( ref = 0 )</td>
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
<td>8 [ 0.05 , 11 ]</td>
<td>8 [ 0.05 , 11 ]</td>
</tr>
<tr class="even">
<td>n.parameters</td>
<td>2</td>
<td>2</td>
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

## Example 2. Survival analysis

The second example is survival analysis (outcome.type=‘SURVIVAL’) to
estimate the effects on the risk of diabetic retinopathy at 8 years of
follow-up, treating macrovascular complications as censoring. 15
covariates and censoring strata are specified in nuisance.model= and
strata=, respectively.

``` r
data(diabetes.complications)
diabetes.complications$d <- (diabetes.complications$epsilon>0)
output <- polyreg(nuisance.model = Event(t,d) ~ age+sex+bmi+hba1c+diabetes_duration
          +drug_oha+drug_insulin+sbp+ldl+hdl+tg+current_smoker+alcohol_drinker+ltpa, 
          exposure = 'fruitq1', strata='strata', data = diabetes.complications,
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
| median.follow.up     | 8 \[ 0.05 , 11 \]                        |
| n.parameters         | 2                                        |
| optimization.message | Function criterion near zero             |

## Example 3. Binomial analysis

Binomial analysis is conducted if outcome.type=‘BINOMIAL’. Outcomes of
observations censored before 8 years of follow-up are now treated as
complication-free.

``` r
output <- polyreg(nuisance.model = d ~ age+sex+bmi+hba1c+diabetes_duration
          +drug_oha+drug_insulin+sbp+ldl+hdl+tg+current_smoker+alcohol_drinker+ltpa, 
          exposure = 'fruitq1', strata='strata', data = diabetes.complications,
          effect.measure1='RR', time.point=8, outcome.type='BINOMIAL')
msummary(output$summary, statistic = c("conf.int"), exponentiate = TRUE)
```

|                      | event 1 (no competing risk)              |
|----------------------|------------------------------------------|
| fruitq1 ( ref = 0 )  | 1.362                                    |
|                      | \[1.152, 1.611\]                         |
| effect.measure       | RR of fruitq1 ( ref = 0 )                |
| n.events             | 362                                      |
| n.events.exposed     | 114 events in 258 exposed observations   |
| n.events.unexposed   | 248 events in 720 unexposed observations |
| median.follow.up     | 0 \[ 0 , 0 \]                            |
| n.parameters         | 2                                        |
| optimization.message | Function criterion near zero             |

## Example 4. Competing risks analysis

The code below specifies direct polytomous regression of both of
competing events (outcome.type=‘COMPETINGRISK’). Initial values of
regression parameters are imported by data.initial.values.

``` r
data.initial.values <- c(-20.01, 1.591, -0.2065, 1.962, 10.05, 4.004, 0.9624, 1.582, 2.414, -0.02140,
                         1.072, 1.304, -0.6444, -1.600, -0.01816, 0.4460, -29.25, 6.586, -1.090, 3.515,
                         4.048, -5.012, 0.05348, 0.1609, 6.432, -0.4783, -3.230, -4.135, 1.451, 0.04964,
                         3.822, -0.08824)
output <- polyreg(nuisance.model = Event(t,epsilon) ~ age+sex+bmi+hba1c+diabetes_duration
          +drug_oha+drug_insulin+sbp+ldl+hdl+tg+current_smoker+alcohol_drinker+ltpa, 
          exposure = 'fruitq1', strata='strata', data=diabetes.complications,
          effect.measure1='RR', time.point=8, outcome.type='COMPETINGRISK',
          data.initial.values=data.initial.values)
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
<td>1.557</td>
<td>0.913</td>
</tr>
<tr class="even">
<td></td>
<td>[1.311, 1.850]</td>
<td>[0.461, 1.805]</td>
</tr>
<tr class="odd">
<td>effect.measure</td>
<td>RR of fruitq1 at 8 ( ref = 0 )</td>
<td>RR of fruitq1 at 8 ( ref = 0 )</td>
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
<td>8 [ 0.05 , 11 ]</td>
<td>8 [ 0.05 , 11 ]</td>
</tr>
<tr class="even">
<td>n.parameters</td>
<td>3</td>
<td>3</td>
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
