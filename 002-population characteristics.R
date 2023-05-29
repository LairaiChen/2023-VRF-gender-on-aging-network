setwd('...')
source('codes/0_final/000-functions_and_paths.R')
path.result <- 'results/demo_by_sex'
pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x) # data with group info
  g <- factor(rep(1:length(x), times=sapply(x, length))) # e.g. 11222
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test or ANOVA
    if (length(unique(g)) == 2) {
      p <- t.test(y ~ g)$p.value
    } else {
      p <- broom::tidy(aov(y ~ g))$p.value[[1]]
    }
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

statistic.and.pval <- function(x, ...) {
  y <- unlist(x) # data with group info
  g <- factor(rep(1:length(x), times=sapply(x, length))) # e.g. 11222
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test or ANOVA
    if (length(unique(g)) == 2) {
      temp <- t.test(y ~ t)
      p <- temp$p.value
      statistic <- temp$statistic %>% unname()
    } else {
      temp <- broom::tidy(aov(y ~ g))
      p <- temp$p.value[[1]]
      statistic <- temp$statistic[[1]]
    }
  } else {
    # For categorical variables, perform a chi-squared test of independence
    temp <- chisq.test(table(y, g))
    p <- temp$p.value
    statistic <- temp$statistic %>% unname()
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  # c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
  #str_c(round(statistic, 3), '\n(', case_when(
  str_c(round(statistic, 3), '(', case_when(
    p < .0001 ~ '****',
    p < .001 ~ '***',
    p < .01 ~ '**',
    p < .05 ~ '*',
    p < .1 ~ '&',
    TRUE ~ 'N.S'
  ), ')')
}


data <- reduce(
  list(
    # age, sex, and education
    read_csv(path.demo),
    read_csv(path.risk.binary) %>% mutate(
      eid = eid,
      `Over drink`= alcohol.bin,
      `Ever smoke` = smoke.bin,
      Obesity = obesity.bin,
      Hypertensive = hypertension.bin,
      Diabetic = diabetes.bin,
      .keep = 'none'
    ),
    read_csv(path.risk.continuous) %>%
      mutate(
        eid = eid,
        `Pack years`=smoke,
        `Alcohol consumption`=alcohol,
        BMI = BMI,
        `Systolic pressure` = systolic,
        `Diastolic pressure` = diastolic,
        `Pulse pressure` = systolic-diastolic, # Ref: EHJ-2019-Simon Cox
        HbA1c = hba1c, # https://labs.selfdecode.com/blog/hba1c-what-is-it-and-why-is-it-important/#:~:text=The%20normal%20range%20for%20HbA1c,%2Fmol)%20may%20indicate%20diabetes.
        #Triglyceride = triglycerides,
        .keep='none'), # Risk factors
    # Brain IDP, global metric and cognition
    read_csv(path.brain.idp) %>% dplyr::mutate(eid, `Brain volume`= BrainV_raw/1000000, GM=GMV_raw/1000000, WM=WMV_raw/1000000, WMH=WMHV_raw/1000, .keep='none'),
    read_csv(path.gene.risk) %>% dplyr::mutate(eid = eid, ApoE=ifelse(`ApoE 4 status` == 'carrier', 1, 0), .keep='none'),
    read_csv(path.global.gretna.wide) %>% dplyr::select(eid, `Global efficiency` = eg, `Local efficiency`=eloc),
    read_csv(path.cognition.raw)
  ), full_join) %>% 
  filter(!is.na(`Global efficiency`))

table1::label(data$age) <- 'Age'
table1::label(data$sex) <- 'Sex'
table1::label(data$education) <- 'Education'

table1::units(demo$age) <- "years"
table1::units(demo$education) <- "years"



tbl1.html <-  table1::table1(~age + education+
                               `Over drink`+`Ever smoke` + Obesity + Hypertensive + Diabetic +
                               `Pack years` + `Alcohol consumption` + BMI + `Systolic pressure` + `Diastolic pressure` + `Pulse pressure` + HbA1c + ApoE +
                               `Brain volume` + GM + WM + WMH +
                               `Global efficiency` + `Local efficiency` + 
                               `fluid intelligence` + `reaction time` + `pairs matching`+ `TMT-A`+ `TMT-B`+ `numeric memory`+ `symbol digit`+ `tower rearrange`+ `matrix pattern` | sex, data = data %>% mutate(sex = as.factor(ifelse(sex == 1, 'Male', 'Female'))) %>% filter(!is.na(sex)), extra.col = list(`Statistics`=statistic.and.pval))


tbl1.df <- as_tibble(tbl1.html)
write_csv(tbl1.df,str_c(path.result, '/demo_2w.csv'))
