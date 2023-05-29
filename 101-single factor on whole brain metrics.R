setwd('/Users/lairaichin/Desktop/Research_Lairai/2022-Sex_difference_of_normal_aging_WM_topology/')
source('codes/000-functions_and_paths.R')
path.result <- 'results/single_factor_on_global_measurements'
create_when_absent(path.result)

# T tests on residuals -----------------------
data <- reduce(
  list(
    # age, sex, and education
    read_csv(path.demo) %>% mutate(education=ifelse(education=='Colleged',1,0)),
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
  ), full_join)

data.long <- data %>% 
  pivot_longer(`Brain volume`:`matrix pattern`, names_to = 'measure', values_to = 'Y') %>% 
  pivot_longer(Alcohol:Diabetes, names_to = 'VRF', values_to = 'group') %>% 
  mutate(group = ifelse(group == 0, 'No', 'Yes')) %>% 
  na.omit()

data.nest <- data.long %>% 
  group_by(measure, VRF, sex) %>% 
  nest() %>% 
  mutate(`Sample size` = map_dbl(data, function(df) {
    nrow(df %>% na.omit() %>% filter(abs(Y - mean(Y)) <= 3 * sd(Y)))
  })) 


result.t <- data.nest %>% 
  mutate(map_dfr(data, function(df) {
    df <- df %>% 
      filter(abs(Y - mean(Y)) <= 3 * sd(Y)) %>% 
      mutate(Y = Y %>% fit_Y_from_X(age) %>% fit_Y_from_X(age ^ 2) %>% fit_Y_from_X(education) %>% fit_Y_from_X(ApoE))
    coef.matrix <- t.test(Y~group, data = df) %>% broom::tidy()
    return(
      tibble(
        p.value = coef.matrix %>% pull(p.value),
        t = coef.matrix %>% pull(statistic),
        df = coef.matrix %>% pull(parameter)
      )
    )
  })) %>% 
  mutate(
    significance = case_when(
      p.value<.0001 ~ '****',
      p.value<.001 ~ '***',
      p.value<.01 ~ '**',
      p.value<.05 ~ '*',
      p.value<.1 ~ '&',
      p.value<=1 ~ 'NS',
      TRUE ~ NA_character_
    )
  ) %>% 
  select(-data)
write_csv(result.t, str_c(path.result, '/T_results.csv'))

# Visualization -----------------------------------------------------------

# https://stackoverflow.com/a/45614547
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

data <- reduce(
  list(
    # age, sex, and education
    read_csv(path.demo) %>% mutate(education=ifelse(education=='Colleged',1,0)),
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
                                               .keep='none'),
    read_csv(path.global.gretna.wide) %>% dplyr::select(eid, `Global efficiency` = eg, `Local efficiency`=eloc),
    read_csv(path.cognition.raw)
  ), full_join)

data.long <- data %>% 
  pivot_longer(c(`Brain volume`:`matrix pattern`, -WMH), names_to = 'measure', values_to = 'Y') %>% 
  pivot_longer(Alcohol:Diabetes, names_to = 'VRF', values_to = 'group') %>% 
  mutate(group = ifelse(group == 0, 'No', 'Yes')) %>% 
  na.omit()


# Fit for age, education and apoe
ggplot(
  data = data.long %>% 
    filter(
      VRF %in% c('Alcohol', 'Smoking', 'Obesity', 'Hypertension', 'Diabetes') &
        measure %in% c('Global efficiency', 'Local efficiency') &
        TRUE) %>% 
    group_by(sex, measure, VRF) %>% 
    filter(abs(Y - mean(Y)) <= 3 * sd(Y)) %>% 
    mutate(Y = Y %>% fit_Y_from_X(age) %>% fit_Y_from_X(age ^ 2) %>% fit_Y_from_X(education) %>% fit_Y_from_X(ApoE)) %>% 
    ungroup() %>% 
    mutate(sex = ifelse(sex == 1, 'Male', 'Female')),
  #mapping = aes(x=group, y=Y, fill=sex)
  mapping = aes(x=sex, y=Y, fill=group)
) +
  geom_split_violin(trim= T,color="white") + 
  geom_boxplot(width = 0.1, 
               alpha = 0.3,
               color = 'white',
               #fill = 'white',
               notch = FALSE, 
               notchwidth = .3, 
               outlier.shape = NA,
               coef=0
               #outlier.size = 0.01
  ) +
  ggsci::scale_fill_futurama(alpha = 0.6) + 
  facet_wrap(measure ~ VRF, scales = 'free') +
  theme_classic() 



