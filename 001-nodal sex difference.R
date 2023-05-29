setwd('...')
source('codes/000-functions_and_paths.R')
path.result <- 'results/sex_difference_nodal'

data <- reduce(
  list(
    read_csv(path.demo) %>% mutate(education=ifelse(education=='Colleged',1,0)),
    read_csv(path.gene.risk) %>% dplyr::mutate(eid = eid, ApoE=ifelse(`ApoE 4 status` == 'carrier', 1, 0), .keep='none'),
    read_csv(path.brain.idp) %>% dplyr::mutate(eid,WMH = WMHV_raw, .keep='none'),
    read_csv(path.nodal.gretna) %>% filter(measure %in% c('NE', 'NLE'))
  ), full_join)


data.long <- data %>% 
  mutate(region = str_sub(region, start = str_length('region') + 1) %>% as.numeric()) %>% 
  filter(is.finite(value) &!is.na(education) & !is.na(sex) & is.finite(ApoE)) %>% 
  na.omit() %>% 
  group_by(sex, measure, region) %>% 
  # mutate(value = value %>% fit_Y_from_X(age) %>% fit_Y_from_X(age^2) %>% fit_Y_from_X(education) %>% fit_Y_from_X(ApoE)) %>% 
  mutate(value = value %>% fit_Y_from_X(age) %>% fit_Y_from_X(age^2) %>% fit_Y_from_X(education) %>% fit_Y_from_X(ApoE) %>% fit_Y_from_X(WMH)) %>% 
  ungroup()


data.nest <- 
  data.long %>% 
  select(sex, region, value, measure) %>% 
  group_by(measure, region) %>% 
  nest()

result <- data.nest %>% 
  mutate(map_dfr(data, function(df) {
    df <- df %>% 
      na.omit() %>% filter(abs(value - mean(value)) <= 3 * sd(value)) %>% 
      mutate(sex = ifelse(sex==1, 'Male', 'Female'))
    coef.matrix <- t.test(value~sex, data = df) %>% broom::tidy()
    return(
      tibble(
        p.value = coef.matrix %>% pull(p.value),
        t = coef.matrix %>% pull(statistic),
        df = coef.matrix %>% pull(parameter)
      )
    )
  })) %>% select(-data) %>% 
  group_by(measure) %>% 
  mutate(p.adj = p.adjust(p.value, method = 'bonferroni')) %>% 
  ungroup() 

  
result <- full_join(result, read_csv(path.description.bna) %>% select(region, `full name`), by = 'region')
write_csv(result, str_c(path.result, '/result.csv'))
# Saves the results for brain net viewer -------------------------------------------------------
result <- read_csv(str_c(path.result, '/result.csv'))

result.adj <- result %>% mutate(t.adj = ifelse(p.adj <= 0.05, abs(t), 0))
bna.description <- read_csv(path.description.bna)
result.nest <- result.adj %>% group_by(measure) %>% nest()

path.result.bnv <- str_c(path.result, '/bnv')
create_when_absent(path.result.bnv)
walk2(result.nest$measure, result.nest$data, function(measure, data) {
  sum_nodes_bnv(
    data,
    bna.description,
    save.size = T,
    size.col = 't.adj',
    path.node.bnv = str_c(path.result.bnv, '/sex_difference_', measure, '.node'),
    path.node.csv = str_c(path.result.bnv, '/sex_difference_', measure, '.csv'),
  )
  data %>% arrange(region) %>% pull('t.adj') %>% 
    write_lines(str_c(path.result.bnv, '/sex_difference_', measure, '.txt'))
})