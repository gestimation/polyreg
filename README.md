
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

The code below conducts unadjusted analysis focusing on cumulative
incidence probabilities of event 1 and 2 at 24 years. Regression
coefficients and variance covariance matrix of both exposure (fruitq1)
and covariates (intercept in this case) in the models are presented.

``` r
library(polyreg)
data(diabetes.complications)
output <- polyreg(nuisance.model = Event(t,epsilon) ~ 1, exposure = 'fruitq1',
          cens.model = Event(t,epsilon==0) ~ 1, data = diabetes.complications,
          effect.measure1='RR', effect.measure2='RR', time.point=8, outcome.type='C')
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
<col style="width: 28%" />
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
<td>0.251</td>
<td>0.018</td>
</tr>
<tr class="even">
<td></td>
<td>[0.212, 0.296]</td>
<td>[0.014, 0.025]</td>
</tr>
<tr class="odd">
<td>fruitq1</td>
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
<td>279 events in 978 observations</td>
<td>79 events in 978 observations</td>
</tr>
<tr class="odd">
<td>median.follow.up</td>
<td>7.9973 [ 0.0493 , 11.0034 ]</td>
<td>7.9973 [ 0.0493 , 11.0034 ]</td>
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
<td>1.50821400358489e-11</td>
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

The second example is survival analysis (outcome.type=‘SURVIVAL’) to
estimate the effects on the risk of diabetic retinopathy at 8 years of
follow-up, treating macrovascular complications as censoring. 15
covariates and censoring strata are specified in nuisance.model and
cens.model, respectively.

``` r
data(diabetes.complications)
diabetes.complications$d <- (diabetes.complications$epsilon>0)
output <- polyreg(nuisance.model = Event(t,d) ~ age+sex+bmi+hba1c
          +diabetes_duration+drug_oha+drug_insulin+sbp+ldl+hdl+tg
          +current_smoker+alcohol_drinker+ltpa, exposure = 'fruitq1',
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
output <- polyreg(nuisance.model = Event(t,epsilon) ~ age+sex+bmi+hba1c
          +diabetes_duration+drug_oha+drug_insulin+sbp+ldl+hdl+tg
          +current_smoker+alcohol_drinker+ltpa, exposure = 'fruitq1',
          cens.model=Event(t,epsilon==0)~strata(strata),data=diabetes.complications,
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
