# ISSUES ------------------------------------------------------------------------------------------------------
# https://forum.posit.co/t/understanding-glmnets-multinomial-logistic-regression-coefficients-not-k-1/52031
# https://stats.stackexchange.com/questions/121155/glmnet-how-to-make-sense-of-multinomial-parameterization
# -------------------------------------------------------------------------------------------------------------

setwd("~/surfdrive/LogitMDA/rcm-association")
library(glmnet)
load("~/surfdrive/LogitMDA/rcm-association/ISSPdata.RData")

# --------------------------------------------------------------------
# --------------------------------------------------------------------
# ISSUE
# 10fold cv is throwing out rows of the contingency table
# thats NOT ok as when
# LT is thrown out only 959 participants are left out
# ZA is thrown out only 3022 participants are left out
# --------------------------------------------------------------------
# --------------------------------------------------------------------

# out = glmnet(x = Z, y = N, family = "multinomial")
# plot(out)
# check stability of 10-fold CV  procedure:
# Amin = A1se = list()
# for(r in 1:10){
#   cv.out <- cv.glmnet(x = Z, y = N, family = "multinomial")
#   Bmin = coef(out, s = cv.out$lambda.min)
#   B1se = coef(out, s = cv.out$lambda.1se)
#   Bmin = cbind(Bmin$G, Bmin$Po, Bmin$NPo, Bmin$Ro, Bmin$Fa)
#   B1se = cbind(B1se$G, B1se$Po, B1se$NPo, B1se$Ro, B1se$Fa)
#   colnames(Bmin) = colnames(B1se) =colnames(N)
#   Amin[[r]] = Bmin  
#   A1se[[r]] = B1se  
# }
# Amin
# A1se

# just once
# cv.out <- cv.glmnet(x = Z, y = N, family = "multinomial", keep = TRUE)
# cv.out$foldid
# plot(cv.out)
# Bmin = coef(out, s = cv.out$lambda.min); B1se = coef(out, s = cv.out$lambda.1se)
# Bmin = cbind(Bmin$G, Bmin$Po, Bmin$NPo, Bmin$Ro, Bmin$Fa); Bmin
# B1se = cbind(B1se$G, B1se$Po, B1se$NPo, B1se$Ro, B1se$Fa); B1se


# --------------------------------------------------------------------
# making LONG matrices
# --------------------------------------------------------------------
YY = matrix(0, 44847, 5)
XX = matrix(NA, 44847, 10)
tel = 0
for(r in 1:34){
  for(c in 1:5){
    nrc = N[r,c]
    if(nrc > 0){
      XX[(tel+1):(tel + nrc), ] = outer(rep(1,nrc), Z[r, -1, drop = TRUE])
      YY[(tel+1):(tel +  nrc), c] = 1
      tel = tel + nrc
    }
  }
}
colnames(YY) = colnames(N)
colnames(XX) = colnames(Z)[-1]

out = glmnet(x = XX, y = YY, family = "multinomial")
plot(out)

cv.out <- cv.glmnet(x = XX, y = YY, family = "multinomial", keep = TRUE)
#cv.out$foldid
plot(cv.out)

cv.out2 <- cv.glmnet(x = XX, y = YY, family = "multinomial", type.measure = "mse", keep = TRUE)
plot(cv.out2)

Bmin = coef(out, s = cv.out$lambda.min); B1se = coef(out, s = cv.out$lambda.1se)
Bmin = cbind(Bmin$G, Bmin$Po, Bmin$NPo, Bmin$Ro, Bmin$Fa); Bmin
B1se = cbind(B1se$G, B1se$Po, B1se$NPo, B1se$Ro, B1se$Fa); B1se
colnames(Bmin) = colnames(B1se) =colnames(N)

xtable::xtable(as.matrix(B1se))
PIhat = predict(out, newx = Z[, -1], s = cv.out$lambda.1se, type = "response")
PIhat = drop(PIhat)
Nhat = diag(rowSums(N)) %*% PIhat
sum(((Nhat - N)^2/Nhat))


ZX = scale(X)
mx = attr(ZX,"scaled:center")
sdx = attr(ZX,"scaled:scale")

# all predictors at mean, LFPf from 0 to 100 %
Xnew.lfpf = outer(rep(1, 11), colMeans(X)); Xnew.lfpf[, 2] = seq(0, 100, by = 10)
Xnew.lfpf = scale(Xnew.lfpf, center = mx, scale = sdx)
Xnew.lfpf # check - correct

round(predict(out, newx = Xnew.lfpf[, -1], s = cv.out$lambda.1se, type = "response"), 3)
#           G    Po   NPo    Ro    Fa
# [1,]  0.392 0.435 0.072 0.012 0.089
# [2,]  0.512 0.319 0.070 0.011 0.087
# [3,]  0.627 0.219 0.063 0.010 0.081
# [4,]  0.725 0.142 0.054 0.008 0.071
# [5,]  0.802 0.088 0.044 0.007 0.059
# [6,]  0.859 0.053 0.035 0.005 0.048
# [7,]  0.900 0.031 0.027 0.004 0.038
# [8,]  0.929 0.018 0.020 0.003 0.029
# [9,]  0.949 0.010 0.015 0.002 0.023
# [10,] 0.964 0.006 0.012 0.002 0.017
# [11,] 0.974 0.003 0.009 0.001 0.013

PIhat = predict(out, newx = Xnew.lfpf[, -1], s = cv.out$lambda.1se, type = "response")
piplot2 = data.frame(probs = matrix(PIhat, ncol = 1), 
                    resp = rep(colnames(PIhat), each = 11),
                    x = rep(seq(0,100, by = 10), 5))


ggplot(data = piplot2, aes(x = x, y = probs, group= resp)) + 
  geom_line(aes(col= resp)) + geom_point(aes(col= resp)) +
  xlab("Labor Force Participation of Females") + 
  ylab("Predicted Profile") + ylim(-0.1,1.1) + 
  scale_x_continuous(breaks=seq(0, 100, 10)) + 
  scale_y_continuous(breaks=seq(0, 1, .1)) + 
  theme_bw()


# From ISSPanalysis.R you can get piplot, which can be combined with piplot2
piplotc = rbind(piplot, piplot2)
piplotc$method = c(rep("nsca", 55), rep("glmnet", 55))

ggplot(data = piplotc, aes(x = x, y = probs, group=interaction(resp, method))) + 
  geom_line(aes(col= resp, linetype = method)) + geom_point(aes(col= resp)) +
  xlab("Labor Force Participation of Females") + 
  ylab("Predicted Profile") + ylim(-0.1,1.1) + 
  scale_x_continuous(breaks=seq(0, 100, 10)) + 
  scale_y_continuous(breaks=seq(0, 1, .1)) + 
  theme_bw()

# predictions for India with varying levels of LaborForceParticipation of Females (observed value = 21.8%)
Xnew16.lfpf = outer(rep(1, 11), X[16, ]); Xnew16.lfpf[, 2] = seq(0, 100, by = 10)
Xnew16.lfpf = scale(Xnew16.lfpf, center = mx, scale = sdx)
Xnew16.lfpf # check - correct

PIhat = predict(out, newx = Xnew16.lfpf, s = cv.out$lambda.1se, type = "response")
#          G    Po   NPo    Ro    Fa
# [1,] 0.368 0.384 0.082 0.020 0.146
# [2,] 0.479 0.280 0.079 0.018 0.144
# [3,] 0.587 0.192 0.072 0.016 0.133
# [4,] 0.683 0.125 0.062 0.013 0.117
# [5,] 0.762 0.078 0.051 0.011 0.098
# [6,] 0.823 0.048 0.040 0.008 0.080
# [7,] 0.870 0.028 0.032 0.006 0.064
# [8,] 0.904 0.016 0.024 0.005 0.050
# [9,] 0.929 0.009 0.018 0.004 0.039
# [10,] 0.948 0.005 0.014 0.003 0.030
# [11,] 0.962 0.003 0.010 0.002 0.023
piplot2 = data.frame(probs = matrix(PIhat, ncol = 1), 
                     resp = rep(colnames(PIhat), each = 11),
                     x = rep(seq(0,100, by = 10), 5))


plt2 = ggplot(data = piplot2, aes(x = x, y = probs, group= resp)) + 
  geom_line(aes(col= resp)) + geom_point(aes(col= resp)) +
  xlab("Labor Force Participation of Females") + 
  ylab("Predicted Profile") + 
  scale_x_continuous(breaks=seq(0, 100, 10)) + 
  scale_y_continuous(breaks=seq(0, 1, .1)) + 
  theme_bw()

# From ISSPanalysis.R you can get piplot, which can be combined with piplot2
piplotc = rbind(piplot, piplot2)
piplotc$method = c(rep("nsca", 55), rep("glmnet", 55))
piplotc$method <- factor(piplotc$method, levels=c("nsca", "glmnet"))

plt3 = ggplot(data = piplotc, aes(x = x, y = probs, group= resp)) + 
  geom_line(aes(col= resp)) + geom_point(aes(col= resp)) +
  facet_grid(cols = vars(method)) + 
  xlab("Labor Force Participation of Females") + 
  ylab("Predicted Profile") + 
  scale_x_continuous(breaks=seq(0, 100, 10)) + 
  scale_y_continuous(breaks=seq(0, 1, .1)) + 
  theme_bw()
plt3

ggsave("~/surfdrive/LogitMDA/rcm-association/figures_issp/predictions.pdf", plot = plt3, width = 11.7, height = 6.0, units = "in", limitsize = FALSE)
