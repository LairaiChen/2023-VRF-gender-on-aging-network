setwd('...')
source('codes/000-functions_and_paths.R')
data <- reduce(list(
  read_csv(path.demo) %>% dplyr::select(eid, age, sex),
  read_csv(path.brain.idp) %>% dplyr::mutate(eid, 
                                             `Brain volume`=BrainV_raw, 
                                             GMV = GMV_raw, 
                                             WMV = WMV_raw,
                                             WMH = WMHV_raw,
                                             # WMH=log(WMHV_raw/BrainV_raw), 
                                             .keep='none'),
  read_csv(path.global.gretna.wide) %>% dplyr::select(eid, `Global integration` = eg, `Global segregation`=eloc)
), full_join) %>% 
  pivot_longer(cols = `Brain volume`:`Global segregation`, names_to = 'phenotype', values_to = 'Y') %>% 
  na.omit()

# Scale the variables
df <- data %>% 
  mutate(sex=ifelse(sex==1, 'Male', 'Female')) %>% mutate(sex=as.factor(sex)) %>% 
  filter(between(age, 50, 80)) %>% 
  #filter(phenotype %in% c('Global integration', 'Global segregation', 'GMV', 'WMV'))
  filter(phenotype %in% c('Global integration', 'Global segregation', 'GMV', 'WMV', 'WMH')) 

ggplot(data= df, 
       aes(x=age, y=Y, color = phenotype)) +
  geom_smooth(aes(linetype = sex), method='loess', alpha = 0.1) +
  ggsci::scale_color_locuszoom() + # 4
  theme_classic()  

# Gets the fitted data
p <- ggplot(data= df, 
            aes(x=age, y=Y, color = phenotype)) +
  geom_smooth(aes(linetype = sex), method='loess', alpha = 0.1) +
  ggsci::scale_color_locuszoom() + # 4
  #ggsci::scale_color_startrek() +
  theme_classic()  



derivative <- function(df) {
  res <- tibble(x = rep(0, nrow(df) - 1), y = rep(0, nrow(df) - 1))
  for(i in 1:length(df$x) - 1) {
    res$x[i] <- df$x[i]
    res$y[i] <- (df$y[i+1] - df$y[i])/(df$x[i+1] - df$x[i])
  }
  return(res)
}

pg <- ggplot_build(p)
df.fit <- pg$data[[1]] %>% dplyr::select(colour, linetype, x, y)
df.res.nest <- df.fit %>% group_by(colour, linetype) %>% nest()
df.res <- df.res.nest %>% mutate(
  map_dfr(data, function(df){
    derive.first <- derivative(df)
    derive.second <- derivative(derive.first)
    derive.first <- mutate(derive.first, y = abs(y))
    derive.second <- mutate(derive.second, y = abs(y))
    return(
      tibble(
        `Quickest decline age` = derive.first$x[which(derive.first$y == max(derive.first$y))],
        `Quickest decline`= max(derive.first$y),
        `Quickest acceleration age` = derive.second$x[which(derive.second$y == max(derive.second$y))],
        `Quickest acceleration` = max(derive.second$y)
      )
    )
  })
) %>% 
  dplyr::select(-data)

write_csv(df.res, 'results/age_trajectories_gam/age_derivatives.csv')