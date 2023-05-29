setwd('...')
source('codes/000-functions_and_paths.R')
path.result <- 'results/combined_VRS_on_nodal' %>% create_when_absent()

data <- reduce(
  list(
    read_csv(path.demo) %>% mutate(education = ifelse(education == 'Colleged', 1, 0)),
    read_csv(path.gene.risk) %>% dplyr::mutate(eid = eid, ApoE=ifelse(`ApoE 4 status` == 'carrier', 1, 0), .keep='none'),
    read_csv(path.risk.binary) %>% mutate(
      eid = eid,
      Alcohol = alcohol.bin,
      Smoking = smoke.bin,
      Obesity = obesity.bin,
      Hypertension = ifelse(hypertension.bin == T, 1, 0),
      Diabetes = ifelse(diabetes.bin == T, 1, 0),
      .keep = 'none'
    ),
    read_csv(path.brain.idp) %>% dplyr::mutate(eid, WMH = WMHV_raw, .keep='none'),
    read_csv(path.nodal.gretna) %>% filter(measure %in% c('NE', 'NLE'))
  ),
  full_join
) %>% 
  mutate(
    VRS = case_when(
      sex == 1 ~ 0.495 * Obesity + 0.25 * Alcohol + 0.355 * Smoking + 0.53 * Diabetes + 0.515 * Hypertension,
      sex == 0 ~ 0.435 * Obesity + 0.065 * Alcohol + 0.165 * Smoking + 0.415 * Diabetes + 0.715 * Hypertension,
      TRUE ~ NA_real_
    )) %>% 
  select(-(Alcohol:Diabetes))

data.long <- data %>% 
  dplyr::rename(Y = value) %>% 
  filter(is.finite(Y)) %>% 
  na.omit() %>% 
  mutate(region = str_sub(region, start = str_length('region') + 1) %>% as.numeric())


# Both genders------------------------------------------------------------

result <- data.long %>% 
  group_by(region, measure) %>% 
  nest() %>% 
  mutate(
    `Sample size` = map_dbl(data, function(df) {
      nrow(df %>% na.omit() %>% filter(abs(Y - mean(Y)) <= 3 * sd(Y)))
    })) %>% 
  mutate(
    map_dfr(data, function(df) {
      df <- df %>% mutate(
        age = scale.lairai(age), 
        VRS = scale.lairai(VRS),
        Y = scale.lairai(Y)
      )
      coef.matrix <-  lm(Y ~ age + I(age^2) + sex + education +  age:sex + ApoE + sex + VRS + sex:VRS, df) %>% broom::tidy()
      return(
        tibble(
          p.value.vrs = coef.matrix %>% filter(term == 'VRS') %>% pull(p.value),
          estimate.vrs = coef.matrix %>% filter(term == 'VRS') %>% pull(estimate),
          se.vrs = coef.matrix %>% filter(term == 'VRS') %>% pull(std.error),
          
          p.value.sex = coef.matrix %>% filter(term == 'sex') %>% pull(p.value),
          estimate.sex = coef.matrix %>% filter(term == 'sex') %>% pull(estimate),
          se.sex = coef.matrix %>% filter(term == 'sex') %>% pull(std.error),
          
          p.value.interaction = coef.matrix %>% filter(term == 'sex:VRS') %>% pull(p.value),
          estimate.interaction = coef.matrix %>% filter(term == 'sex:VRS') %>% pull(estimate),
          se.interaction = coef.matrix %>% filter(term == 'sex:VRS') %>% pull(std.error)
        )
      )
    })
  ) %>% select(-data)

result.adj <- result %>% 
  group_by(measure) %>% 
  mutate(
    p.vrs.fdr = p.adjust(p.value.vrs, method = 'fdr'),
    p.sex.fdr = p.adjust(p.value.sex, method = 'fdr'),
    p.interaction.fdr = p.adjust(p.value.interaction, method = 'fdr')
  ) %>% 
  ungroup() %>% 
  full_join(read_csv(path.description.bna) %>% select(region, `full name`), by = 'region')

write_csv(result.adj, str_c(path.result, '/results_both_gender.csv'))


# Single sex--------------------------------------------------------------
result <- data.long %>% 
  group_by(sex, region, measure) %>% 
  nest() %>% 
  mutate(
    `Sample size` = map_dbl(data, function(df) {
      nrow(df %>% na.omit() %>% filter(abs(Y - mean(Y)) <= 3 * sd(Y)))
    })) %>% 
  mutate(
    map_dfr(data, function(df) {
      #df <- data.nest$data[[1]]
      df <- df %>% mutate(
        age = scale.lairai(age), 
        VRS = scale.lairai(VRS),
        Y = scale.lairai(Y)
      )
      coef.matrix <-  lm(Y ~ age + I(age^2) + education + ApoE + VRS, df) %>% broom::tidy()
      return(
        tibble(
          p.value.vrs = coef.matrix %>% filter(term == 'VRS') %>% pull(p.value),
          estimate.vrs = coef.matrix %>% filter(term == 'VRS') %>% pull(estimate),
          se.vrs = coef.matrix %>% filter(term == 'VRS') %>% pull(std.error)
        )
      )
    })
  ) %>% select(-data)

result.adj <- result %>% 
  group_by(sex, measure) %>% 
  mutate(
    p.vrs.fdr = p.adjust(p.value.vrs, method = 'fdr'),
    p.sex.fdr = p.adjust(p.value.sex, method = 'fdr'),
    p.interaction.fdr = p.adjust(p.value.interaction, method = 'fdr')
  ) %>% 
  ungroup() %>% 
  full_join(read_csv(path.description.bna) %>% select(region, `full name`), by = 'region')

write_csv(result, str_c(path.result, '/results_single_gender.csv'))

