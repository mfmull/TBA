############################################################
# Border-proximate groundwater irrigation analysis
#
# Purpose:
#   - Classify transboundary aquifer dyads by irrigation patterns
#   - Compare current vs potential future classifications
#   - Produce figures and underlying datasets
#
# Inputs:
#   FinalDiads.csv
#
# Outputs:
#   ./plts/*.pdf figures
#   ./plts/*.csv datasets underlying figures
#   SI2.csv classification table
############################################################


############################################################
# 1. Libraries
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(multcompView)
library(purrr)
library(tibble)
library(ggalluvial)
library(vcd)
library(grid)



############################################################
# 2. Parameters
############################################################

eps    <- 1e-3
thresh <- 0.99

cls <- rev(c("IU0","IU1","UB","DA","BL"))

cols <- c(
  IU0 = "#D3D3D3",
  IU1 = "#A6CEE3",
  UB  = "#6A3D9A",
  DA  = "#FF7F00",
  BL  = "#E31A1C"
)



############################################################
# 3. Load data
############################################################

diads <- read.csv("FinalDiads.csv") %>%
  distinct(code, CC_1, CC_2, .keep_all = TRUE)



############################################################
# 4. Current dyad classification
############################################################

diadsP <- diads %>%
  mutate(
    rG1_1 = GW_B1_1 / GW_1,
    rG1_2 = GW_B1_2 / GW_2,
    rG2_1 = GW_B22_1 / GW_1,
    rG2_2 = GW_B22_2 / GW_2,
    
    class = case_when(
      
      pmax(GW_B22_1, GW_B22_2) < eps &
        pmax(UR_B22_1, UR_B22_2) < eps ~ "IU0",
      
      (GW_B22_1 < eps & UR_B22_1 < eps) |
        (GW_B22_2 < eps & UR_B22_2 < eps) ~ "IU1",
      
      pmin(GW_B22_1, GW_B22_2) < eps ~ "UB",
      
      pmax(rG1_1, rG1_2) > thresh ~ "BL",
      
      pmax(rG2_1, rG2_2) > thresh ~ "DA",
      
      TRUE ~ "DA"
    ),
    
    class = factor(class, levels = cls),
    GW = GW_1 + GW_2
  )



############################################################
# 5. Future classification
############################################################

diadsF <- diads %>%
  mutate(
    GW_1 = pmax(GW_1, CR3_1),
    GW_2 = pmax(GW_2, CR3_2),
    GW_B22_1 = pmax(GW_B22_1, CR3_B22_1),
    GW_B22_2 = pmax(GW_B22_2, CR3_B22_2),
    GW_B1_1 = pmax(GW_B1_1, CR3_B1_1),
    GW_B1_2 = pmax(GW_B1_2, CR3_B1_2)
  ) %>%
  mutate(
    rG1_1 = GW_B1_1 / GW_1,
    rG1_2 = GW_B1_2 / GW_2,
    rG2_1 = GW_B22_1 / GW_1,
    rG2_2 = GW_B22_2 / GW_2,
    
    classF = case_when(
      
      pmax(GW_B22_1, GW_B22_2) < eps &
        pmax(UR_B22_1, UR_B22_2) < eps ~ "IU0",
      
      (GW_B22_1 < eps & UR_B22_1 < eps) |
        (GW_B22_2 < eps & UR_B22_2 < eps) ~ "IU1",
      
      pmin(GW_B22_1, GW_B22_2) < eps ~ "UB",
      
      pmax(rG1_1, rG1_2) > thresh ~ "BL",
      
      pmax(rG2_1, rG2_2) > thresh ~ "DA",
      
      TRUE ~ "DA"
    ),
    
    classF = factor(classF, levels = cls)
  )



############################################################
# 6. Export classification table
############################################################

diadsC <- diadsP %>%
  select(code, name, CC_1, CC_2, class) %>%
  left_join(
    diadsF %>% select(code, CC_1, CC_2, classF),
    by = c("code","CC_1","CC_2")
  )




############################################################
# 7. Alluvial transition plot
############################################################

p1 <- ggplot(diadsC, aes(axis1 = class, axis2 = classF)) +
  geom_alluvium(aes(fill = class), width = 1/12) +
  geom_stratum(width = 1/6, aes(fill = class), color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Current","Future"), expand = c(.05,.05)) +
  scale_fill_manual(values = cols) +
  theme_minimal()

ggsave("./plts/B200/cls_1.pdf",p1,width=4,height=3)


data_alluvial <- diadsC %>%
  count(class,classF,name="n")

write.csv(data_alluvial,
          "./plts/B200/data_cls_1_alluvial.csv",
          row.names=FALSE)



############################################################
# 8. Tukey tests and faceted bar plots
############################################################

diadsP <- diadsP %>%
  mutate(
    FIR = (GW_1 + GW_2)/(GW_C_1 + GW_C_2),
    FBW = (GW3_1 + GW3_2)/(GW_1 + GW_2),
    GW  = (GW_1 + GW_2)*1000
  )

d_long <- diadsP %>%
  filter(GW>0) %>%
  select(class,FIR,FBW) %>%
  pivot_longer(-class,names_to="Variable",values_to="Value")


tukey_df <- d_long %>%
  group_by(Variable) %>%
  nest() %>%
  mutate(res = map(data,\(d){
    
    fit <- aov(Value~class,d)
    tuk <- TukeyHSD(fit,conf.level=.9)
    
    letters <- multcompLetters4(fit,tuk,threshold=.1)$class[[1]] |>
      tibble::enframe("class","letters")
    
    d %>%
      group_by(class) %>%
      summarise(
        mean = mean(Value,na.rm=TRUE),
        se   = sd(Value,na.rm=TRUE)/sqrt(n()),
        ci   = list(t.test(Value)$conf.int),
        .groups="drop"
      ) %>%
      mutate(
        ci_low  = map_dbl(ci,1),
        ci_high = map_dbl(ci,2)
      ) %>%
      select(-ci) %>%
      left_join(letters)
    
  })) %>%
  select(-data) %>%
  unnest(res) %>%
  mutate(
    Variable=factor(Variable,c("FGW","FBW")),
    class=factor(class,levels=rev(cls))
  )


ptuck <- ggplot(tukey_df,aes(class,mean,fill=class))+
  geom_col(width=.6)+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=.1)+
  geom_text(aes(label=letters,
                y=ci_high+0.1*(ci_high-ci_low)),size=3.5)+
  geom_hline(yintercept=0,linetype="dashed")+
  facet_wrap(~Variable,scales="free_y",nrow=1)+
  scale_fill_manual(values=cols,name="")+
  labs(x=NULL,y="Mean value")+
  theme_classic()+
  theme(
    legend.position="top",
    strip.background=element_blank(),
    strip.text=element_text(face="bold")
  )

ggsave("./plts/B200/cls_2.pdf",ptuck,width=5,height=3)


data_tukey <- tukey_df %>%
  select(Variable,class,mean,se,letters)

write.csv(data_tukey,
          "./plts/B200/data_cls_2_bars.csv",
          row.names=FALSE)


