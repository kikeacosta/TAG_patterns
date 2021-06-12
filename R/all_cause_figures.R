
library(tidyverse)

IN <- read.table("Output/EstDiffRateAging.txt",sep=" ", header = TRUE) %>% 
  as.tibble()

# geom_ribbon(aes(ymin = lower, ymax = upper), color = "#0055AA40")

# 1) all-cause example:
# subset a country baseline and 2020 eta (log mort)
# show both with confidence bands.

colnames(IN)
DE <- 
  IN %>% 
  filter(country == "Germany") %>% 
  mutate(eta.1.point = exp(eta.1),
         eta.1.low = exp(eta.1.low),
         eta.1.up = exp(eta.1.up),
         eta.2.point = exp(eta.2),
         eta.2.low = exp(eta.2.low),
         eta.2.up = exp(eta.2.up)) %>% 
  select(age, eta.1.point,eta.1.low,eta.1.up,
         eta.2.point,eta.2.low,eta.2.up) %>% 
  pivot_longer(eta.1.point:eta.2.up, values_to = "rate") %>% 
  separate(name, sep = "\\.", into = c(NA,"year","which")) %>% 
  pivot_wider(names_from = "which", values_from = rate) %>% 
  mutate(year = recode(year, `1` = "2020", `2` = "baseline"))
library(scales)
DE %>% 
  ggplot(aes(x=age, y = point, ymin = low, ymax = up, color = year))+
  geom_line() +
  geom_ribbon(fill = "#00555560") + 
  scale_y_log10(labels = comma) +
  guides(color = FALSE) + 
  labs(y = "log mortality rate", x = "Age", title = "Germany (total)") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        title=element_text(size=14)
        ) + 
  annotate("text",x=50,y=.01,label="2020",size=12)+ 
  annotate("text",x=75,y=.003,label="baseline",size=12)

# GR <- 
#   IN %>% 
#   filter(country == "Liechtenstein") %>% 
#   mutate(eta.1.point = exp(eta.1),
#          eta.1.low = exp(eta.1.low),
#          eta.1.up = exp(eta.1.up),
#          eta.2.point = exp(eta.2),
#          eta.2.low = exp(eta.2.low),
#          eta.2.up = exp(eta.2.up)) %>% 
#   select(age, eta.1.point,eta.1.low,eta.1.up,
#          eta.2.point,eta.2.low,eta.2.up) %>% 
#   pivot_longer(eta.1.point:eta.2.up, values_to = "rate") %>% 
#   separate(name, sep = "\\.", into = c(NA,"year","which")) %>% 
#   pivot_wider(names_from = "which", values_from = rate) %>% 
#   mutate(year = recode(year, `1` = "2020", `2` = "baseline"))
# 
# GR %>% 
#   ggplot(aes(x=age, y = point, ymin = low, ymax = up, color = year))+
#   geom_line() +
#   geom_ribbon(fill = "#00555560") + 
#   scale_y_log10(labels = comma) +
#   guides(color = FALSE) + 
#   labs(y = "log mortality rate", x = "Age", title = "Greece (total)") +
#   theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=14),
#         title=element_text(size=14)
#   ) + 
#   annotate("text",x=50,y=.01,label="2020",size=12)+ 
#   annotate("text",x=75,y=.003,label="baseline",size=12)



# 2) rate of aging 
# same country, baseline and 2020 with CIs

DE <- 
  IN %>% 
  filter(country == "Germany") %>% 
  select(age, roa.1.point = roa.1, roa.1.low,roa.1.up,
         roa.2.point = roa.2,roa.2.low,roa.2.up) %>% 
  pivot_longer(roa.1.point:roa.2.up, values_to = "rate") %>% 
  separate(name, sep = "\\.", into = c(NA,"year","which")) %>% 
  pivot_wider(names_from = "which", values_from = rate) %>% 
  mutate(year = recode(year, `1` = "2020", `2` = "baseline"))

DE %>% 
  ggplot(aes(x=age, y = point, ymin = low, ymax = up, color = year, fill = year))+
  geom_line() +
  geom_ribbon(alpha=.5) + 
  scale_y_log10(labels = comma) +
  guides(color = FALSE, fill = FALSE) + 
  labs(y = "rate of aging", x = "Age", title = "Germany (total)") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        title=element_text(size=14)
  ) + 
  annotate("text",x=25,y=.1,label="baseline",size=12)+ 
  annotate("text",x=86,y=.03,label="2020",size=12)

# 3) difference in rate of aging 2020-baseline, with CIs

IN %>% 
  filter(country == "Germany") %>% 
  ggplot(aes(x=age, y = delta, ymin = delta.low, ymax = delta.up))+
  geom_line() +
  geom_ribbon(alpha=.5, color = "#00555560") + 
  labs(y = "difference in rate of aging", x = "Age", title = "Germany (total)") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        title=element_text(size=14)) + 
  geom_hline(yintercept = 0, color = "#FF000080")
  # ) + 

  # annotate("text",x=25,y=.1,label="baseline",size=12)+ 
  # annotate("text",x=86,y=.03,label="2020",size=12)

# 4) plot nr (3) facetted for all countries, make ribbon color vary depending on whether bounds include 0.



pad_sig_ribbons <- function(chunk){
  chunk %>% 
    mutate(check = sig != lag(sig)) %>% 
    dplyr::filter(check) %>% 
    select(-check) %>% 
    mutate(sig = !sig) %>% 
    bind_rows(chunk) %>% 
    arrange(age)
}

IN %>% 
  #group_by(country) %>% 
  # mutate(sig = delta.low >= 0 | delta.up <= 0) %>% 
  # do(pad_sig_ribbons(chunk = .data)) %>% 
  # ungroup() %>% 
  #dplyr::filter(country %in% c("Spain","Germany")) %>% 
  ggplot(aes(x=age,y=delta,ymax=delta.low,ymin=delta.up))+
  geom_line() +
  geom_ribbon(fill = "#00555560") + 
  geom_hline(yintercept = 0, color = "#FF000080") +
  #scale_fill_manual(values = c("#00555540","#DD220070")) + 
  facet_wrap(~country) +
  ylim(-.15,.15) + 
  coord_cartesian(clip = "off") +
  labs(y = "difference in rate of aging", x = "Age") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        title=element_text(size=14)) +
  geom_rect(data = dplyr::filter(IN, country %in% 
                                   c("Moldova","Poland","Ukraine","Armenia","Georgia","Romania","Montenegro","Greece","Croatia","Estonia","Bulgaria","Finland","France")), 
            fill = NA, colour = "red", xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf) 
 

IN %>% 
  dplyr::filter(country %in% c("Moldova","Poland","Ukraine","Armenia","Georgia")) %>% 
  #group_by(country) %>% 
  # mutate(sig = delta.low >= 0 | delta.up <= 0) %>% 
  # do(pad_sig_ribbons(chunk = .data)) %>% 
  # ungroup() %>% 
  #dplyr::filter(country %in% c("Spain","Germany")) %>% 
  ggplot(aes(x=age,y=delta,ymax=delta.low,ymin=delta.up))+
  geom_line() +
  geom_ribbon(fill = "#00555560") + 
  geom_hline(yintercept = 0, color = "#FF000080") +
  #scale_fill_manual(values = c("#00555540","#DD220070")) + 
  facet_wrap(~country) +
  ylim(-.25,.25) + 
  coord_cartesian(clip = "off")


