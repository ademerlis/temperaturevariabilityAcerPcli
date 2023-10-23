There are so many ways I can plot the photosynthetic efficiency of the corals in the temperature treatment and the CBASS. I'm not sure which visualization is best, so I'm going to put all of them here for now.

## 1. Raw Fv/Fm during the 28-d temperature treatment period (repeated measures).

<img width="528" alt="Screen Shot 2023-10-22 at 12 42 58 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/9b07e4f5-1883-40e4-833d-41e5a627bc9d">

Stats for raw Fv/Fm:
Things to consider: 
- repeated measures over time (include time as a random or fixed effect?)
- what is the statistical question: is the Fv/Fm decline over time significantly different between treatments?
    - this could be best tested then by converting fv/fm to a rate: fv/fm decline per day.
    - Compare normalized versus unnormalized to see which fits linear model or statistical test best

Genotype: are there significant differences of genotype within species as well? Are some genotypes more influenced by the treatment than others?

<img width="623" alt="Screen Shot 2023-10-22 at 12 48 15 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/3090f7f3-9785-4add-a490-4cbe0ca3c600">


## 2. Decline of Fv/Fm over treatment period, normalized to initial measurement

Plot points as scatterplot and then do line of best fit and get slope of line and test significance of slope 

<img width="630" alt="Screen Shot 2023-08-28 at 12 32 19 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/d1036d46-b024-4a49-84e2-590b24975de0">

## 3. Spread of normalized Fv/Fm data at the end of the treatment period.

<img width="623" alt="Screen Shot 2023-10-13 at 2 08 30 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/3e20fbc9-5639-4568-a594-bbb21f667f2b">


## 4. Reduction of Fv/Fm over 28-d Treatment, normalized to initial measurement (all genotypes combined)

<img width="621" alt="Screen Shot 2023-10-22 at 3 08 57 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/ff574741-5dfd-45dc-a2cc-0c5a015fd521">


Stats for this:

**Acer stats**

```{r}
## Model comparison function
#two *nested* models as input, in this case obj1 is the model without the predictor variable and obj2 is the model with the predictor variable 
lrt <- function (obj1, obj2) {
  L0 <- logLik(obj1)
  L1 <- logLik(obj2)
  L01 <- as.vector(- 2 * (L0 - L1))
  df <- attr(L1, "df") - attr(L0, "df")
  list(L01 = L01, df = df,
       "p-value" = pchisq(L01, df, lower.tail = FALSE))
}

ipam_acer_raw %>% drop_na(fvfm_loss) -> ipam_acer_raw

#model 1 = full model Treatment + (1|Colony) + (1|Tank)
glmm_Acer_full <- glmmTMB(fvfm_loss_norm ~ Treatment + (1|Colony) + (1|Tank), family=beta_family(link = "logit"), data=ipam_acer_raw)

summary(glmm_Acer_full) 

qqnorm(resid(glmm_Acer_full))
qqline(resid(glmm_Acer_full))

boxplot(resid(glmm_Acer_full)~ipam_acer_raw$Treatment)
boxplot(resid(glmm_Acer_full)~ipam_acer_raw$Colony)
boxplot(resid(glmm_Acer_full)~ipam_acer_raw$Tank)

#model 2 = Treatment + (1|Colony)
Acer_fvfm_GLMM_Colony <- glmmTMB::glmmTMB(fvfm_loss_norm ~ Treatment + (1|Colony), family=beta_family(link = "logit"), data=ipam_acer_raw)
summary(Acer_fvfm_GLMM_Colony)

qqnorm(resid(Acer_fvfm_GLMM_Colony)) 
qqline(resid(Acer_fvfm_GLMM_Colony))

boxplot(resid(Acer_fvfm_GLMM_Colony)~ipam_acer_raw$Treatment) 
boxplot(resid(Acer_fvfm_GLMM_Colony)~ipam_acer_raw$Colony)

#model 3 = Treatment + (1|Tank)
Acer_fvfm_GLMM_Tank <- glmmTMB::glmmTMB(fvfm_loss_norm ~ Treatment + (1|Tank), family=beta_family(link = "logit"), data=ipam_acer_raw)
summary(Acer_fvfm_GLMM_Tank)

qqnorm(resid(Acer_fvfm_GLMM_Tank)) 
qqline(resid(Acer_fvfm_GLMM_Tank))

boxplot(resid(Acer_fvfm_GLMM_Tank)~ipam_acer_raw$Treatment) 
boxplot(resid(Acer_fvfm_GLMM_Tank)~ipam_acer_raw$Colony)

lrt(Acer_fvfm_GLMM_Colony, glmm_Acer_full) #significant

lrt(Acer_fvfm_GLMM_Tank, glmm_Acer_full) #significant

# GLMM Acer Full: AIC =   -366.2
# GLMM Acer Colony: AIC =  -335.9
# GLMM Acer Tank: AIC =  -358.2

#based on AIC and LRT p-values (lowest AIC is best, if p-value is significant that means there is a significant difference in the full model versus the reduced model), the full model is the best model.

Anova(glmm_Acer_full, type = "II")

emmeans(glmm_Acer_full, pairwise ~ Treatment)

summary(glht(glmm_Acer_full, linfct=mcp(Treatment="Tukey")), test = adjusted(type = "bonferroni"))
```

<img width="446" alt="Screen Shot 2023-10-22 at 2 59 40 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/a42662ab-5e79-4889-b173-9d173c49cb47">

<img width="455" alt="Screen Shot 2023-10-22 at 3 00 22 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/94af28bf-b24f-45e4-be30-1b2eaac5c655">

<img width="511" alt="Screen Shot 2023-10-22 at 3 00 01 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/93f4ded5-7d9a-42df-93be-597cf74afbcf">


**Pcli stats**

```{r}
ipam_pcli_raw %>% drop_na(fvfm_loss) -> ipam_pcli_raw

ipam_pcli_raw %>% 
  filter(fvfm_loss_norm > 0) -> ipam_pcli_raw

#model 1 = full model Treatment*time_point + (1|Colony) + (1|Tank)
glmm_pcli_full <- glmmTMB(fvfm_loss_norm ~ Treatment + (1|Colony) + (1|Tank), family=beta_family(link = "logit"), data=ipam_pcli_raw)

summary(glmm_pcli_full) 

qqnorm(resid(glmm_pcli_full))
qqline(resid(glmm_pcli_full))

boxplot(resid(glmm_pcli_full)~ipam_pcli_raw$Treatment)
boxplot(resid(glmm_pcli_full)~ipam_pcli_raw$Colony)
boxplot(resid(glmm_pcli_full)~ipam_pcli_raw$Tank)

#model 2 = Treatment*time_point + (1|Colony)
pcli_fvfm_GLMM_Colony <- glmmTMB::glmmTMB(fvfm_loss_norm ~ Treatment + (1|Colony), family=beta_family(link = "logit"), data=ipam_pcli_raw)
summary(pcli_fvfm_GLMM_Colony)

qqnorm(resid(pcli_fvfm_GLMM_Colony)) 
qqline(resid(pcli_fvfm_GLMM_Colony))

boxplot(resid(pcli_fvfm_GLMM_Colony)~ipam_pcli_raw$Treatment) 
boxplot(resid(pcli_fvfm_GLMM_Colony)~ipam_pcli_raw$Colony)

#model 3 = Treatment*Date + (1|Tank)
pcli_fvfm_GLMM_Tank <- glmmTMB::glmmTMB(fvfm_loss_norm ~ Treatment + (1|Tank), family=beta_family(link = "logit"), data=ipam_pcli_raw)
summary(pcli_fvfm_GLMM_Tank)

qqnorm(resid(pcli_fvfm_GLMM_Tank)) 
qqline(resid(pcli_fvfm_GLMM_Tank))

boxplot(resid(pcli_fvfm_GLMM_Tank)~ipam_pcli_raw$Treatment) 
boxplot(resid(pcli_fvfm_GLMM_Tank)~ipam_pcli_raw$Colony)

lrt(pcli_fvfm_GLMM_Colony, glmm_pcli_full) #significant

lrt(pcli_fvfm_GLMM_Tank, glmm_pcli_full) #significant

# GLMM Full: AIC =   -466.4 
# GLMM Colony: AIC = -374.0
# GLMM Tank: AIC =  -458.9

#based on AIC and LRT p-values (lowest AIC is best, if p-value is significant that means there is a significant difference in the full model versus the reduced model), the full model is the best model.

Anova(glmm_pcli_full, type = "II")

emmeans(glmm_pcli_full, pairwise ~ Treatment)

summary(glht(glmm_pcli_full, linfct=mcp(Treatment="Tukey")), test = adjusted(type = "bonferroni"))

```

<img width="385" alt="Screen Shot 2023-10-23 at 11 21 54 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/c03c972c-3b9e-41e2-ab94-595be6baba4d">

<img width="445" alt="Screen Shot 2023-10-23 at 11 22 11 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/20e74ed9-d6d4-4e8c-8305-6426fea4d768">

<img width="515" alt="Screen Shot 2023-10-23 at 11 22 26 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/ce41fb56-263a-4561-8c4d-58870eee76c9">



## 5. Reduction of Fv/Fm over 28-d Treatment, normalized to initial measurement (separating genotypes)

<img width="628" alt="Screen Shot 2023-10-13 at 3 02 47 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/ae8f543c-ef95-47bf-bb7b-a2c91634cf5e">

Stats for this:

**Acer stats**

```{r}
#Creating interactions for post hoc comparisons
ipam_acer_raw$colony_treatment<-interaction(ipam_acer_raw$Colony, ipam_acer_raw$Treatment)
ipam_pcli_raw$colony_treatment<-interaction(ipam_pcli_raw$Colony, ipam_pcli_raw$Treatment)

## Model comparison function

#two *nested* models as input, in this case obj1 is the model without the predictor variable and obj2 is the model with the predictor variable 
lrt <- function (obj1, obj2) {
  L0 <- logLik(obj1)
  L1 <- logLik(obj2)
  L01 <- as.vector(- 2 * (L0 - L1))
  df <- attr(L1, "df") - attr(L0, "df")
  list(L01 = L01, df = df,
       "p-value" = pchisq(L01, df, lower.tail = FALSE))
}

#Acer GLM with Colony as a co-variate

#model 1 = full model colony_treatment + (1|Tank)
glmm_Acer_full <- glmmTMB(fvfm_loss_norm ~ colony_treatment + (1|Tank), family=beta_family(link = "logit"), data=ipam_acer_raw)

summary(glmm_Acer_full) 

qqnorm(resid(glmm_Acer_full))
qqline(resid(glmm_Acer_full))

boxplot(resid(glmm_Acer_full)~ipam_acer_raw$Treatment)
boxplot(resid(glmm_Acer_full)~ipam_acer_raw$Colony)
boxplot(resid(glmm_Acer_full)~ipam_acer_raw$Tank)

#model 2 = Treatment*time_point + (1|Colony)
Acer_fvfm_GLMM_notank <- glmmTMB::glmmTMB(fvfm_loss_norm ~ colony_treatment, family=beta_family(link = "logit"), data=ipam_acer_raw)
summary(Acer_fvfm_GLMM_notank)

qqnorm(resid(Acer_fvfm_GLMM_notank)) 
qqline(resid(Acer_fvfm_GLMM_notank))

boxplot(resid(Acer_fvfm_GLMM_notank)~ipam_acer_raw$Treatment) 
boxplot(resid(Acer_fvfm_GLMM_notank)~ipam_acer_raw$Colony)

lrt(Acer_fvfm_GLMM_notank, glmm_Acer_full) #significant

# GLMM Acer Full: AIC =   -371.8
# GLMM Acer no tank: AIC =  -342.3

#based on AIC and LRT p-values (lowest AIC is best, if p-value is significant that means there is a significant difference in the full model versus the reduced model), the full model is the best model.

Anova(glmm_Acer_full, type = "II")

emmeans(glmm_Acer_full, pairwise ~ colony_treatment)

summary(glht(glmm_Acer_full, linfct=mcp(Treatment="Tukey")), test = adjusted(type = "bonferroni"))
```

<img width="628" alt="Screen Shot 2023-10-22 at 3 11 37 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/177501ba-1a5b-4e90-8f8c-6972b8ced055">


<img width="549" alt="Screen Shot 2023-10-22 at 3 11 09 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/65205f4e-f91f-4a73-9534-acaef4e64a96">

<img width="513" alt="Screen Shot 2023-10-22 at 3 11 23 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/dbf558ef-8ca8-49df-a452-4b74bcd14710">


**Pcli stats**

```{r}
#model 1 = full model colony_treatment + (1|Tank)
glmm_pcli_full <- glmmTMB(fvfm_loss_norm ~ colony_treatment + (1|Tank), family=beta_family(link = "logit"), data=ipam_pcli_raw)

summary(glmm_pcli_full) 

qqnorm(resid(glmm_pcli_full))
qqline(resid(glmm_pcli_full))

boxplot(resid(glmm_pcli_full)~ipam_pcli_raw$Treatment)
boxplot(resid(glmm_pcli_full)~ipam_pcli_raw$Colony)
boxplot(resid(glmm_pcli_full)~ipam_pcli_raw$Tank)

#model 2 = Treatment*time_point + (1|Colony)
pcli_fvfm_GLMM_notank <- glmmTMB::glmmTMB(fvfm_loss_norm ~ colony_treatment, family=beta_family(link = "logit"), data=ipam_pcli_raw)
summary(pcli_fvfm_GLMM_notank)

qqnorm(resid(pcli_fvfm_GLMM_notank)) 
qqline(resid(pcli_fvfm_GLMM_notank))

boxplot(resid(pcli_fvfm_GLMM_notank)~ipam_pcli_raw$Treatment) 
boxplot(resid(pcli_fvfm_GLMM_notank)~ipam_pcli_raw$Colony)

lrt(pcli_fvfm_GLMM_notank, glmm_pcli_full) #significant

# GLMM Acer Full: AIC =   -468.0 
# GLMM Acer no tank: AIC =  -375.1

#based on AIC and LRT p-values (lowest AIC is best, if p-value is significant that means there is a significant difference in the full model versus the reduced model), the full model is the best model.

Anova(glmm_pcli_full, type = "II")

emmeans(glmm_pcli_full, pairwise ~ colony_treatment)

summary(glht(glmm_pcli_full, linfct=mcp(colony_treatment="Tukey")), test = adjusted(type = "bonferroni"))
```

<img width="425" alt="Screen Shot 2023-10-23 at 11 25 46 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/b2c8d9b3-f27f-41f9-aecc-fe7728360593">

<img width="487" alt="Screen Shot 2023-10-23 at 11 26 15 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/03fcb0bc-6b12-4625-8e61-9916f57c8a9b">

<img width="476" alt="Screen Shot 2023-10-23 at 11 26 33 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/5784b042-5502-4f40-89a0-3dd15ac18aba">


## 6. Loss of fv/fm per day normalized to initial measurement

<img width="619" alt="Screen Shot 2023-10-22 at 1 34 03 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/bf71964f-0da6-4ea1-afd2-aa715b2aa34e">

## 7. Loss of fv/fm per day (not normalized to initial measurement)

<img width="616" alt="Screen Shot 2023-10-22 at 1 36 08 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/7962059f-a782-4246-a76e-32633b760614">

Numbers 5-7 all look the same. 

## 7. Dose-Response Curves Following Acute Heat Stress using CBASS

1. Normalized to initial Fv/Fm 
<img width="624" alt="Screen Shot 2023-08-28 at 5 05 43 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/2d0ffe76-d3e0-4ce2-b0b5-32c154d39706">


