setwd('...')
source('codes/000-functions_and_paths.R')
path.result <- create_when_absent('results/single_factor_on_nodal')

path.save <- create_when_absent(str_c(path.result, '/T_and_barplots'))

data <- reduce(
  list(
    # age, sex, and education
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
    read_csv(path.nodal.gretna) %>% filter(measure %in% c('NE', 'NLE'))
  ), full_join
)

data.long <- data %>% 
  pivot_longer(Alcohol:Diabetes, names_to = 'VRF', values_to = 'X') %>% 
  filter(is.finite(value) & !is.na(X) & !is.na(education) & !is.na(sex)) %>% 
  na.omit()

data.fit <-  data.long %>% 
  group_by(measure, sex, VRF, region) %>% 
  filter(abs(value - mean(value, na.rm = T)) <= 3*sd(value, na.rm = T)) %>% 
  mutate(value = value %>% fit_Y_from_X(age)%>% fit_Y_from_X(age^2) %>% fit_Y_from_X(education) %>% fit_Y_from_X(ApoE)) %>% 
  ungroup() %>% 
  mutate(
    sex = factor(sex, levels = c(0, 1), labels = c('Female', 'Male')),
    X = factor(X, levels = c(0, 1), labels = c('Not risk', 'At risk'))
  )

df <- data.fit %>% 
  mutate(group = case_when(
    sex == 'Male' & X == 'Not risk' ~ 'Male without risk',
    sex == 'Male' & X == 'At risk' ~ 'Male at risk',
    sex == 'Female' & X == 'At risk' ~ 'Female at risk',
    sex == 'Female' & X == 'Not risk' ~ 'Female without risk'
  )) %>% 
  mutate(
    group = factor(group, levels = c('Male without risk', 'Male at risk', 'Female without risk', 'Female at risk'), ordered = T)
  ) %>% select(region, value, measure, VRF, group)

df.nest <- df %>% group_by(measure, VRF, region) %>% nest()

comparisons <- list(c('Male without risk', 'Female without risk'), 
                    c('Male at risk', 'Female at risk'),
                    c('Male without risk', 'Male at risk'),
                    c('Female without risk', 'Female at risk'))

bna.label <-  read_csv(path.description.bna) %>% pull('full name')

signif_bonferroni <- function(p, n.compare = 4) {
  p <- p * n.compare
  res <- case_when(
    p < 0.0001 ~ 'p < 0.0001',
    p < 0.001 ~ 'p < 0.001',
    p < 0.01 ~ 'p < 0.01',
    p < 0.05 ~ 'p < 0.05',
    p < 0.01 ~ 'p < 0.1',
    TRUE ~ 'Not significant'
  )
  return (res)
}

future::plan(future::multisession, workers = 5)
furrr::future_pwalk(
  list(df.nest$measure, df.nest$VRF, df.nest$region, df.nest$data),
  function(measure, VRF, region, df) {
    create_when_absent(str_c(path.save, '/', measure))
    create_when_absent(str_c(path.save, '/', measure, '/', VRF))
    p <-
      ggpubr::ggboxplot(df %>% na.omit(), x = 'group', y = 'value', fill = 'group', color = 'black',
                        width = 0.6, # numeric value between 0 and 1 specifying box width.
                        xlab = '',
                        ylab = bna.label[region %>% str_sub(start = 7) %>% as.integer()],
                        outlier.shape = 21,
                        size = 0.2, # Numeric value (e.g.: size = 1). change the size of points and outlines.
                        #font.label = list(size = 10), # Not working
                        palette = c("#FC7E8B","#D96C77","#91A0F2","#828ED9")) +
      ggsignif::geom_signif(comparisons = list(c('Male without risk', 'Female without risk'), 
                                               c('Male at risk', 'Female at risk'),
                                               c('Male without risk', 'Male at risk'),
                                               c('Female without risk', 'Female at risk')),
                            test = 't.test',
                            #map_signif_level = function(p) sprintf("p = %.2g", p),
                            map_signif_level = signif_bonferroni,
                            size = 0.2, # line size
                            textsize = 3,
                            tip_length = 0,
                            step_increase = 0.06,
                            color = 'black'
      ) +
      theme(
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)
      ) +
      guides(fill = 'none')
    
    ggsave(str_c(path.save, '/', measure, '/', VRF, '/', region, '.pdf'), p, width = 4, height = 4)
  }
)
