setwd('...')
source('codes/000-functions_and_paths.R')
path.result <- create_when_absent('results/single_factor_on_modular_measurements')

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
    read_csv(path.modular.gretna)
  ),
  full_join
) 
data.long <- data %>% 
  pivot_longer(Alcohol:Diabetes, names_to = 'VRF', values_to = 'X') %>% 
  pivot_longer(eg:eloc, names_to = 'topology', values_to = 'Y') %>% 
  na.omit()

data.nest <- data.long %>% 
  group_by(module, VRF, topology, sex) %>% 
  nest() %>% 
  mutate(`Sample size` = map_dbl(data, function(df) {
    nrow(df %>% na.omit() %>% filter(abs(Y - mean(Y)) <= 3 * sd(Y)))
  })) 

result.t <- data.nest %>% 
  mutate(map_dfr(data, function(df) {
    df <- df %>% 
      filter(abs(Y - mean(Y)) <= 3 * sd(Y)) %>% 
      mutate(Y = Y %>% fit_Y_from_X(age) %>% fit_Y_from_X(age ^ 2) %>% fit_Y_from_X(education) %>% fit_Y_from_X(ApoE))
    coef.matrix <- t.test(Y~X, data = df) %>% broom::tidy()
    return(
      tibble(
        p.value = coef.matrix %>% pull(p.value),
        t = coef.matrix %>% pull(statistic),
        df = coef.matrix %>% pull(parameter)
      )
    )
  })) %>%
  select(-data) %>% 
  group_by(VRF, topology, sex) %>% 
  mutate(p.bonf = p.adjust(p.value, method = 'bonferroni')) %>% 
  ungroup() %>% 
  mutate(
    signifcance = case_when(
      p.bonf <.0001 ~ '****',
      p.bonf <.001 ~ '***',
      p.bonf <.01 ~ '**',
      p.bonf <.05 ~ '*',
      p.bonf <.1 ~ '&',
      p.bonf <=1 ~ 'NS',
      TRUE ~ NA_character_
    )
  )

write_csv(result.t, str_c(path.result, '/T_results.csv'))

# Lollipop visualization
data.result <- 
  read_csv(str_c(path.result, '/T_results.csv')) %>% 
  mutate(sex = ifelse(sex == 1, 'Male', 'Female')) %>% 
  mutate(significance = factor(signifcance, c('****', '***', '**', '*', '&', 'NS'))) %>% 
  select(sex, VRF, module, topology, t, significance) 

data.lollipop <- data.result %>% 
  select(sex, VRF, module, topology, t) %>% 
  pivot_wider(names_from = topology, values_from = t) %>% 
  mutate(module = factor(module,  levels = c("LM", "VN","SUB", "FPN",  "VAN", "DMN", "DAN", "SMN")))

segment_size <- 0.8
point_size <- 1
grid_size <- 0.1

walk2(c('Male', 'Female'), c('dashed', 'solid'), function(gender, linetype){
  plt <- ggplot() +
    ggsci::scale_color_lancet() +
    geom_segment(
      data = data.lollipop %>% filter(sex == gender),
      aes(x=module, xend=module, y=eloc, yend=eg, color = module), 
      linetype = linetype, size = segment_size) +
    
    ggnewscale::new_scale_color() +
    scale_shape_manual(values=c(8, 19, 15, 17, 18, 4))+
    
    geom_point(
      data = data.result %>% filter(sex == gender),
      aes(x=module, y=t, fill = topology, color = topology),
      shape = 21,
      #fill = 'white',
      stroke = point_size / 3,
      size = point_size * 2
    ) +  
    scale_color_manual(values = c('#619BC0', '#62C17B')) +
    scale_fill_manual(values = c('#619BC0', '#62C17B')) +
    geom_point(
      data = data.result %>% filter(sex == gender),
      aes(x=module, y=t, shape = significance),
      color = 'white',
      stroke = point_size / 3,
      size = point_size
    ) +
    coord_flip()+
    #theme_classic() +
    ylim(c(-3.5, 6.5)) +
    xlab("") +
    ylab("T") +
    theme(
      panel.grid.minor=element_line(colour='lightgray', size = grid_size),
      panel.grid.major=element_line(colour='lightgray', size = grid_size),
      panel.grid.major.x = element_blank(),
      panel.background=element_rect(fill='white')) +
    facet_wrap(~VRF, ncol = 1)
  
  ggsave(
    str_c('results/graphs_modular/single_factor_', gender, '.pdf'),
    plt,
    width = 5, height = 11
  )
})