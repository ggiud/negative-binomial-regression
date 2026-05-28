# APPLICAZIONE TESI NEGATIVE BINOMIAL GIUDITTA ADEZIO

library(tidyverse)
setwd("C:\\Users\\giudi\\OneDrive\\Desktop\\uni\\tesi")

dati <- read.csv("covid data.csv")

dati <- select(dati, - c(note, source, country_code, continent))    

dati2 <- dati %>% 
  group_by(year_week, country) %>%
  pivot_wider(names_from = indicator, values_from = c(weekly_count, rate_14_day, 
                                            cumulative_count), names_sep = '_')
dati2 <- as.data.frame(dati2)

library(visdat)
vis_dat(dati2) +
  scale_fill_manual(values= c('blue', 'seagreen2', 'violet', 'darkorchid'))

dati2 <- dati2[complete.cases(dati2),]

library(ggcorrplot)
ggcorrplot(cor(dati2[,-c(1,3)]), lab = TRUE,                 
           lab_size = 3, colors = c("gold1", "seagreen2", "blue2"),  
           ggtheme = theme_minimal()) 


d1 <- ggplot(dati2, aes(x = weekly_count_deaths)) +
  geom_histogram(aes(y = ..density..), binwidth = 10, fill = "blue2", alpha = 0.7, color = 'white') +
  geom_density(aes(y = ..density..), color = "seagreen2", fill =  "seagreen2", alpha = 0.2, linewidth = 1) +
  labs(title = "Distribuzione weekly count deaths dati completi",
       x = "Weekly Count Deaths",
       y = "Density") +
  scale_x_continuous(limits = c(0, 600)) + 
  theme_bw()

dati2$country <- factor(dati2$country)
levels(dati2$country)

summary(dati2)

data_no_eu <- dati2[dati2$country != "EU/EEA (total)", ]

tab <- xtabs(~ weekly_count_deaths, data_no_eu)
ks.test(tab, 'ppois', mean(data_no_eu$weekly_count_deaths))

d2 <- ggplot(data_no_eu, aes(x = weekly_count_deaths)) +
  geom_histogram(aes(y = ..density..), binwidth = 10, fill = "blue2", alpha = 0.7, color = 'white') +
  geom_density(aes(y = ..density..),  color = "seagreen2", fill =  "seagreen2", alpha = 0.2, size = 1) +
  labs(title = "Distribuzione weekly count deaths dati no eu",
       x = "Weekly Count Deaths",
       y = "Density") +
  scale_x_continuous(limits = c(0, 600)) + 
  theme_bw()

library(gridExtra)
grid.arrange(d1, d2)

min_max_normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

dati2_normalized <- as.data.frame(lapply(dati2[, -c(1, 3, 10, 11)], min_max_normalize))
dati2_long <- reshape2::melt(dati2_normalized)

ggplot(dati2_long, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(show.legend = F, color = 'black', linewidth = 0.2, outlier.alpha = 0.5) +
  scale_fill_manual(values = c("seagreen2", "blue", "cyan", "violet", "blue2", "gold1", "darkorchid")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
  labs(title = "Boxplot delle variabili normalizzate", x = "", y = "Valore Normalizzato")


################################################
# M O D E L L I
################################################

# divisione training test in base al tempo
dati2$year <- as.numeric(sub("-.*", "", dati2$year_week))
dati2$week <- as.numeric(sub(".*-", "", dati2$year_week))

# Ordina i dati in base alla colonna 'date_order'
dati2 <- dati2 %>%
  arrange(year, week)

split_index <- floor(0.8 * nrow(dati2))
train_data <- dati2[1:split_index, ]
test_data <- dati2[(split_index + 1):nrow(dati2), ]


# POISSON
mod_poisson <- glm(weekly_count_deaths ~ population + cumulative_count_cases + 
                     weekly_count_cases + cumulative_count_deaths + country + rate_14_day_cases +
                     rate_14_day_deaths, train_data, family = 'poisson')
summary(mod_poisson)

par(mfrow = c(2,2))
plot(mod_poisson)
par(mfrow = c(1,1))

sum(residuals(mod_poisson, type = 'pearson')^2)/mod_poisson$df.residual

library(car)
outlier_test <- outlierTest(mod_poisson)

outlier_indices <- which(outlier_test$bonf.p < 0.05)
dati_no_outlier <- train_data[-outlier_indices, ]

mod_poisson_no_out <- glm(weekly_count_deaths ~ population + cumulative_count_cases + 
                            weekly_count_cases + cumulative_count_deaths + country + rate_14_day_cases +
                            rate_14_day_deaths, dati_no_outlier, family = 'poisson')
summary(mod_poisson_no_out)

par(mfrow = c(2,2))
plot(mod_poisson_no_out)
par(mfrow = c(1,1))

sum(residuals(mod_poisson_no_out, type = 'pearson')^2)/mod_poisson_no_out$df.residual



# NEGATIVE BINOMIAL
library(MASS)
mod_nb <- glm.nb(weekly_count_deaths ~ population + cumulative_count_cases + 
                   weekly_count_cases + cumulative_count_deaths + country + 
                   rate_14_day_cases + rate_14_day_deaths, train_data)

summary(mod_nb)

(phi.hat_nb <- sum(residuals(mod_nb, type = 'pearson')^2)/mod_nb$df.residual)

par(mfrow = c(2,2))
plot(mod_nb)
par(mfrow = c(1,1))

outlierTest(mod_nb)


# QUASIPOISSON
mod_qp <- glm(weekly_count_deaths ~ population + cumulative_count_cases + 
                weekly_count_cases + cumulative_count_deaths + country + rate_14_day_cases +
                rate_14_day_deaths, train_data, family = 'quasipoisson')
summary(mod_qp)

par(mfrow = c(2,2))
plot(mod_qp)
par(mfrow = c(1,1))


# filtriamo il training e il test togliendo dal dataset l'Europa
train_data2 <- filter(train_data, country!= 'EU/EEA (total)')
train_data2$country <- droplevels(train_data2$country)
table(train_data2$country)

test_data2 <- filter(test_data, country!= 'EU/EEA (total)')


# POISSON 2
mod_poi2 <- glm(weekly_count_deaths ~ cumulative_count_cases + 
                  weekly_count_cases + cumulative_count_deaths + country + rate_14_day_cases +
                  rate_14_day_deaths, data = train_data2, family = poisson)
summary(mod_poi2)

par(mfrow = c(2,2))
plot(mod_poi2)
par(mfrow = c(1,1))

sum(residuals(mod_poi2, type = 'pearson')^2)/mod_poi2$df.residual

library(AER)
dispersiontest(mod_poisson)
dispersiontest(mod_poi2)

## plot diagnostici con ggplot a confronto ##
library(broom)

model_data <- augment(mod_poisson)

# 1. Residuals vs Fitted
p1 <- ggplot(model_data, aes(.fitted, .resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = F, color = "seagreen2") +
  labs(title = "Residuals vs Fitted",
       x = "Fitted values",
       y = "Residuals") +
  theme_bw()

# 2. Normal Q-Q
p2 <- ggplot(model_data, aes(sample = .std.resid)) +
  stat_qq(alpha = 0.5) +
  stat_qq_line(linewidth = 1, color = "gold1") +
  labs(title = "Normal Q-Q",
       x = "Theoretical Quantiles",
       y = "Standardized Residuals") +
  theme_bw()

# 3. Scale-Location (Spread-Location)
p3 <- ggplot(model_data, aes(.fitted, sqrt(abs(.std.resid)))) +
  geom_point(alpha = 0.5) +
  geom_smooth( se = F, color = "blue2") +
  labs(title = "Scale-Location",
       x = "Fitted values",
       y = expression(sqrt(abs("Standardized Residuals")))) +
  theme_bw()

# 4. Residuals vs Leverage
p4 <- ggplot(model_data, aes(.hat, .std.resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth( se = F, color = "violet") +
  labs(title = "Res. vs Leverage",
       x = "Leverage",
       y = "Standardized Residuals") +
  theme_bw()

grid.arrange(p1, p2, p3, p4)

model2_data <- augment(mod_poisson_no_out)

# 1. Residuals vs Fitted
p1_2 <- ggplot(model2_data, aes(.fitted, .resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = F, color = "seagreen2") +
  labs(title = "Residuals vs Fitted",
       x = "Fitted values",
       y = "Residuals") +
  theme_bw()

# 2. Normal Q-Q
p2_2 <- ggplot(model2_data, aes(sample = .std.resid)) +
  stat_qq(alpha = 0.5) +
  stat_qq_line(linewidth = 1, color = "gold1") +
  labs(title = "Normal Q-Q",
       x = "Theoretical Quantiles",
       y = "Standardized Residuals") +
  theme_bw()

# 3. Scale-Location (Spread-Location)
p3_2 <- ggplot(model2_data, aes(.fitted, sqrt(abs(.std.resid)))) +
  geom_point(alpha = 0.5) +
  geom_smooth( se = F, color = "blue2") +
  labs(title = "Scale-Location",
       x = "Fitted values",
       y = expression(sqrt(abs("Standardized Residuals")))) +
  theme_bw()

# 4. Residuals vs Leverage
p4_2 <- ggplot(model2_data, aes(.hat, .std.resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth( se = F, color = "violet") +
  labs(title = "Res. vs Leverage",
       x = "Leverage",
       y = "Standardized Residuals") +
  theme_bw()


grid.arrange(p1_2, p2_2, p3_2, p4_2, nrow = 2)
grid.arrange(p1, p2, p3, p4, p1_2, p2_2, p3_2, p4_2, nrow = 2)


# NEGATIVE BINOMIAL 2
mod_nb2 <- glm.nb(weekly_count_deaths ~ cumulative_count_cases + 
                    weekly_count_cases + cumulative_count_deaths + country + rate_14_day_cases +
                    rate_14_day_deaths, train_data2)

summary(mod_nb2)

sum(residuals(mod_nb2, type = 'pearson')^2)/mod_nb2$df.residual

par(mfrow = c(2,2))
plot(mod_nb2)
par(mfrow = c(1,1))

outlierTest(mod_nb2)


# QUASIPOISSON 2
mod_qp2 <- glm(weekly_count_deaths ~ cumulative_count_cases + 
                 weekly_count_cases + cumulative_count_deaths + country + rate_14_day_cases +
                 rate_14_day_deaths, train_data2, family = 'quasipoisson')
summary(mod_qp2)

(phi.hat_qp2 <- sum(residuals(mod_qp2, type = 'pearson')^2)/mod_qp2$df.residual)


par(mfrow = c(2,2))
plot(mod_qp2)
par(mfrow = c(1,1))



# ZINB
library(pscl)
mod_zinb <- zeroinfl(weekly_count_deaths ~ population + cumulative_count_cases + 
                       weekly_count_cases + cumulative_count_deaths + country + rate_14_day_cases +
                       rate_14_day_deaths, data = train_data, dist = "negbin")
summary(mod_zinb)


# ZIP
mod_zip <- zeroinfl(weekly_count_deaths ~ population + cumulative_count_cases + 
                      weekly_count_cases + cumulative_count_deaths + country + rate_14_day_cases +
                      rate_14_day_deaths, data = train_data, dist = "pois")
summary(mod_zip)


# HURDLE
mod_hurdle <- hurdle(weekly_count_deaths ~ population + cumulative_count_cases + 
                       weekly_count_cases + cumulative_count_deaths + country + rate_14_day_cases +
                       rate_14_day_deaths, data = train_data, dist = "negbin")
summary(mod_hurdle)


# Standardizzo le variabili così da vedere se i problemi di multicollinearità vengono risolti
train_data3 <- train_data
train_data3$weekly_count_cases <- scale(train_data3$weekly_count_cases)
train_data3$cumulative_count_deaths <- scale(train_data3$cumulative_count_deaths)
train_data3$rate_14_day_cases <- scale(train_data3$rate_14_day_cases)
train_data3$rate_14_day_deaths <- scale(train_data3$rate_14_day_deaths)
train_data3$population <- scale(train_data3$population)


# ZINB 2
mod_zinb2 <- zeroinfl(weekly_count_deaths ~ population + weekly_count_cases + cumulative_count_deaths + 
                        rate_14_day_cases + rate_14_day_deaths, data = train_data3, dist = "negbin")
summary(mod_zinb2)

# ZIP 2
mod_zip2 <- zeroinfl(weekly_count_deaths ~ population + weekly_count_cases + cumulative_count_deaths + 
                       rate_14_day_cases + rate_14_day_deaths, train_data3, dist = "pois")
summary(mod_zip2)

# HURDLE 2
mod_hurdle2 <- hurdle(weekly_count_deaths ~ population + cumulative_count_cases + 
                        weekly_count_cases + cumulative_count_deaths + rate_14_day_cases +
                        rate_14_day_deaths, data = train_data3, dist = "negbin")
summary(mod_hurdle2)


# AIC e BIC a confronto
tab_confronto <- cbind(c(AIC(mod_poisson), AIC(mod_poisson_no_out), AIC(mod_poi2), AIC(mod_nb), AIC(mod_nb2), 
                         AIC(mod_zip2), AIC(mod_zinb2), AIC(mod_qp2), AIC(mod_hurdle2)),
                       c(BIC(mod_poisson), BIC(mod_poisson_no_out), BIC(mod_poi2), BIC(mod_nb), BIC(mod_nb2),
                         BIC(mod_zip2), BIC(mod_zinb2), BIC(mod_qp2), BIC(mod_hurdle2)))



colnames(tab_confronto) <- c('AIC', 'BIC')
rownames(tab_confronto) <- c('Poisson', 'Poisson no outliers', 'Poisson no UE',
                             'Negative Binomial', 'Negative Binomial no UE',
                             'Zero Inflated Poisson', 'Zero Inflated NegBin',
                             'Quasi-Poisson', 'Hurdle model')
tab_confronto



# GRAFICO DI CONFRONTO

predicted <- predict(mod_nb, newdata = test_data, type = "response")
uno <- ggplot(test_data, aes(x = weekly_count_deaths, y = predicted)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "gold1", linewidth = 1, linetype = 'dashed') +
  labs(title = "Negative Binomial",
       x = "Valori reali (test data)",
       y = "Valori previsti") +
  scale_y_continuous(limits = c(0,1500))+
  scale_x_continuous(limits = c(0,1500))+
  theme_bw()

predicted_poi <- predict(mod_poisson, newdata = test_data, type = "response")
due <- ggplot(test_data, aes(x = weekly_count_deaths, y = predicted_poi)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "blue2", linewidth = 1, linetype = "dashed") +
  labs(title = "Poisson",
       x = "Valori reali (test data)",
       y = "Valori previsti") +
  scale_y_continuous(limits = c(0,1500))+
  scale_x_continuous(limits = c(0,1500))+
  theme_bw()



predicted_qp <- predict(mod_qp, newdata = test_data, type = "response")
tre <- ggplot(test_data, aes(x = weekly_count_deaths, y = predicted_qp)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "seagreen2", linewidth = 1, linetype = "dashed") +
  labs(title = "Quasi-Poisson",
       x = "Valori reali (test data)",
       y = "Valori previsti") +
  scale_y_continuous(limits = c(0,1500))+
  scale_x_continuous(limits = c(0,1500))+
  theme_bw()


predicted_zinb <- predict(mod_zinb, newdata = test_data, type = "response")
quattro <- ggplot(test_data, aes(x = weekly_count_deaths, y = predicted_zinb)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "violet", linewidth = 1, linetype = "dashed") +
  labs(title = "ZINB",
       x = "Valori reali (test data)",
       y = "Valori previsti") +
  scale_y_continuous(limits = c(0,1500))+
  scale_x_continuous(limits = c(0,1500))+
  theme_bw()


predicted_zip <- predict(mod_zip, newdata = test_data, type = "response")
cinque <- ggplot(test_data, aes(x = weekly_count_deaths, y = predicted_zip)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "darkcyan", linewidth = 1, linetype = "dashed") +
  labs(title = "ZIP",
       x = "Valori reali (test data)",
       y = "Valori previsti") +
  scale_y_continuous(limits = c(0,1500))+
  scale_x_continuous(limits = c(0,1500))+
  theme_bw()


predicted_hurdle <- predict(mod_hurdle, newdata = test_data, type = "response")
sei <- ggplot(test_data, aes(x = weekly_count_deaths, y = predicted_hurdle)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "darkorchid", linewidth = 1, linetype = "dashed") +
  labs(title = "Hurdle",
       x = "Valori reali (test data)",
       y = "Valori previsti") +
  scale_y_continuous(limits = c(0,1500))+
  scale_x_continuous(limits = c(0,1500))+
  theme_bw()

grid.arrange(uno, due, tre, quattro, cinque, sei, ncol = 3)


# metriche di confronto
library(Metrics)

mse_poisson  <- mse(test_data$weekly_count_deaths, predicted_poi)
mse_nb       <- mse(test_data$weekly_count_deaths, predicted)
mse_qp       <- mse(test_data$weekly_count_deaths, predicted_qp)
mse_zip      <- mse(test_data$weekly_count_deaths, predicted_zip)
mse_zinb     <- mse(test_data$weekly_count_deaths, predicted_zinb)
mse_hurdle   <- mse(test_data$weekly_count_deaths, predicted_hurdle)

rmse_poisson <- round(rmse(test_data$weekly_count_deaths, predicted_poi),3)
rmse_nb      <- round(rmse(test_data$weekly_count_deaths, predicted),3)
rmse_qp      <- round(rmse(test_data$weekly_count_deaths, predicted_qp),3)
rmse_zip     <- round(rmse(test_data$weekly_count_deaths, predicted_zip),3)
rmse_zinb    <- round(rmse(test_data$weekly_count_deaths, predicted_zinb),3)
rmse_hurdle  <- round(rmse(test_data$weekly_count_deaths, predicted_hurdle),3)

mae_poisson  <- mae(test_data$weekly_count_deaths, predicted_poi)
mae_nb       <- mae(test_data$weekly_count_deaths, predicted)
mae_qp       <- mae(test_data$weekly_count_deaths, predicted_qp)
mae_zip      <- mae(test_data$weekly_count_deaths, predicted_zip)
mae_zinb     <- mae(test_data$weekly_count_deaths, predicted_zinb)
mae_hurdle   <- mae(test_data$weekly_count_deaths, predicted_hurdle)

r2_poisson <- round(cor(test_data$weekly_count_deaths, predicted_poi)^2,3)
r2_nb      <- round(cor(test_data$weekly_count_deaths, predicted)^2,3)
r2_qp      <- round(cor(test_data$weekly_count_deaths, predicted_qp)^2,3)
r2_zip     <- round(cor(test_data$weekly_count_deaths, predicted_zip)^2,3)
r2_zinb    <- round(cor(test_data$weekly_count_deaths, predicted_zinb)^2,3)
r2_hurdle  <- round(cor(test_data$weekly_count_deaths, predicted_hurdle),3)


tab_metriche <- cbind(c(mse_nb, mse_poisson, mse_qp, mse_zinb, mse_zip, mse_hurdle),
                      c(rmse_nb, rmse_poisson, rmse_qp, rmse_zinb, rmse_zip, rmse_hurdle),
                      c(round(mae_nb, 3), round(mae_poisson, 3), round(mae_qp, 3), round(mae_zinb, 3), round(mae_zip, 3), round(mse_hurdle, 3)),
                      c(r2_nb, r2_poisson, r2_qp, r2_zinb, r2_zip, r2_hurdle))                

colnames(tab_metriche) <- c('MSE', 'RMSE', 'MAE', expression(("R")^2))
rownames(tab_metriche) <- c('Negative Binomial', 'Poisson', 'Quasi-Poisson', 
                            'Zero Inflated Negative Binomial','Zero Inflated Poisson',
                            'Hurdle')

tab_metriche

