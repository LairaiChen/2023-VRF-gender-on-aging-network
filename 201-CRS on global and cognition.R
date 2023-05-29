setwd('...')
source('codes/000-functions_and_paths.R')
path.result <- create_when_absent('results/combined_VRS_on_global_and_cognition')

data <- 
  reduce(
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
      read_csv(path.brain.idp) %>% dplyr::mutate(eid, 
                                                 `Brain volume`=BrainV_raw, 
                                                 GMV = GMV_raw, 
                                                 WMV = WMV_raw,
                                                 WMH = WMHV_raw,
                                                 # WMH=log(WMHV_raw/BrainV_raw), 
                                                 .keep='none'),
      read_csv(path.global.gretna.wide) %>% dplyr::select(eid, `Global efficiency` = eg, `Local efficiency`=eloc),
      read_csv(path.cognition.raw)
    ),
    full_join
  ) %>% 
  mutate(
    VRS = case_when(
      sex == 1 ~ 0.495 * Obesity + 0.25 * Alcohol + 0.355 * Smoking + 0.53 * Diabetes + 0.515 * Hypertension,
      sex == 0 ~ 0.435 * Obesity + 0.065 * Alcohol + 0.165 * Smoking + 0.415 * Diabetes + 0.715 * Hypertension,
      TRUE ~ NA_real_
    )) %>% 
  filter(!(eid %in% id_exlude)) %>% 
  select(-(Alcohol:Diabetes))

data.long <- data %>% 
  pivot_longer(cols=c(`Brain volume`:`matrix pattern`, -WMH), names_to = 'measure', values_to = 'Y') %>% 
  filter(is.finite(Y) & !is.na(VRS) & !is.na(education) & !is.na(sex) & is.finite(ApoE)) %>% 
  na.omit()

# Combined models of both genders -----------------------------------------
data.nest <- data.long %>% 
  group_by(measure) %>% 
  nest() %>% 
  mutate(`Sample size` = map_dbl(data, function(df) {
    nrow(df %>% na.omit() %>% filter(abs(Y - mean(Y)) <= 3 * sd(Y)))
  })) 

result <- data.nest %>% 
  mutate(
    map_dfr(data, function(df) {
      #df <- data.nest$data[[1]]
      df <- df %>% mutate(
        age = scale.lairai(age), 
        VRS = scale.lairai(VRS),
        Y = scale.lairai(Y)
      )
      coef.matrix <-  lm(Y ~ age + I(age^2) + sex + education +  age:sex + ApoE + sex + VRS + sex:VRS, data = df) %>% broom::tidy()
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
  ) %>% 
  select(-data)

write_csv(result, str_c(path.result, '/results_both_gender.csv'))



# Separate models by males and females ------------------------------------
result <- data.long %>% 
  group_by(sex, measure) %>% 
  nest() %>% 
  mutate(`Sample size` = map_dbl(data, function(df) {
    nrow(df %>% na.omit() %>% filter(abs(Y - mean(Y)) <= 3 * sd(Y)))
  })) %>% 
  mutate(
    map_dfr(data, function(df) {
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
  ) %>% 
  select(-data)

write_csv(result, str_c(path.result, '/results_single_gender.csv'))






