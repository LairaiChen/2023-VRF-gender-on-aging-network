setwd('...')
source('codes/000-functions_and_paths.R')
path.result <- 'results/combined_VRS_on_modules'
create_when_absent(path.result)

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
    read_csv(path.modular.gretna)
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
  pivot_longer(eg:eloc, names_to = 'topology', values_to = 'Y') %>% 
  filter(is.finite(Y)) %>% 
  na.omit()

# Both genders ------------------------------------------------------------

result <- data.long %>% 
  group_by(topology, module) %>% 
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

write_csv(result, str_c(path.result, '/results_both_gender.csv'))


# Single gender----------------------------------------
result <- data.long %>% 
  group_by(sex, topology, module) %>% 
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

write_csv(result, str_c(path.result, '/results_single_gender.csv'))

# Visualization for each gender group -------------------------------------

path.save <- 'results/graphs_modular/'
result <- 
  readr::read_csv(str_c(path.result, '/results_single_gender.csv')) %>% 
  mutate(sex = ifelse(sex == 1, 'male', 'female'))

circuit_bar <- function(df, measure, ymin=0, ymax=0.1, by=0.02) {

  data <-  df %>% 
    filter(topology == measure) %>%  
    select(module, value = estimate.vrs) %>% 
    mutate(
      value = abs(value),
      module = factor(module, levels = c("LM", "VN","VAN", "DAN", "FPN", "DMN", "SUB", "SMN"))
    ) 
  plt <-
    ggplot(data) +
    # Make custom panel grid
    geom_hline(
      aes(yintercept = y), 
      data.frame(y=seq(from = ymin, to = ymax, by = by)),
      color = "lightgrey"
    ) + 
    geom_col(
      aes(
        x = module,
        y = value,
        fill = module
      ),
      position = "dodge2",
      show.legend = TRUE,
      alpha = .9
    ) +
    
    # Add dots to represent the mean gain
    geom_point(
      aes(
        x = module,
        y = value
      ),
      size = 3,
      color = "gray12"
    ) +
    
    # Lollipop shaft 
    geom_segment(
      aes(
        x = module,
        y = ymin,
        xend = module,
        #ystart = 0,
        yend = ymax
      ),
      linetype = "dashed",
      color = "gray12"
    ) + 
    # Make it circular!
    coord_polar() +
  ggsci::scale_fill_lancet() +
    theme(
      # Remove axis ticks and text
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      # Move the legend to the bottom
      legend.position = "bottom",
      # Make the background white and remove extra grid lines
      panel.background = element_rect(fill = "white", color = "white"),
      panel.grid = element_blank(),
      panel.grid.major.x = element_blank()
    ) 
  return(plt)
}
walk(c('male', 'female'), function(gender){
  walk(c('eg', 'eloc'), function(measure){
    result.data <- result %>% filter(sex == gender)
    ggsave(str_c(path.save, '/combined_', measure, '_', gender, '.pdf'), 
           circuit_bar(result.data, measure),
           width = 4, height = 4)
  })
})
