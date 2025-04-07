#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# load data
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------

# load("ISSPdata.RData")
load("~/surfdrive/LogitMDA/rcm-association/ISSPdata.RData")
Xo = X

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# analysis of full sample
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------

library(CAvariants)
source("~/surfdrive/LogitMDA/rcm-association/nsca.R")

#----------------------------------------------------------------------------
# usual NSCA
#----------------------------------------------------------------------------
out = CAvariants(t(N), catype = "NSCA")
summary(out)
# 
# Total inertia  0.019 
# 
# Inertias, percent inertias and cumulative percent inertias of the row and column space
# 
# inertia inertiapc cuminertiapc
# value1   0.017    86.261       86.261
# value2   0.002    11.714       97.976
# value3   0.000     1.850       99.826
# value4   0.000     0.174      100.000
# 
# Predictability Index for Variants of Non symmetrical Correspondence Analysis:
#   
#   Tau Index predicting from column 
# 
# [1] 0.074
# 
# C-statistic 13239.45 and p-value 0 

output = nsca(N)
cumsum(output$d^2)/sum(output$d^2) * 100
# 86.26117  97.97559  99.82558 100.00000 100.00000
output$loss
# 0.0003886595
sum(output$d[3:5]^2)
# 0.0003886595

#----------------------------------------------------------------------------
# Constrained NSCA
#----------------------------------------------------------------------------
output2 = nsca.it(N, X = Z)
output2$loss
# 0.007404161 

(sum(out$inertias[, 1]) - output2$loss)/ sum(out$inertias[, 1])
# 0.6143395
round(output2$B, digits = 3)

# 3D
output3 = nsca.it(N, X = Z, S = 3)
output3$loss
# 0.007317413
(sum(out$inertias[, 1]) - output3$loss)/ sum(out$inertias[, 1])
# 0.6188579

# 4D
output4 = nsca.it(N, X = Z, S = 4)
output4$loss
# 0.007309155 
(sum(out$inertias[, 1]) - output4$loss)/ sum(out$inertias[, 1])
# 0.619288

#----------------------------------------------------------------------------
# Constrained NSCA with Lasso penalty
#----------------------------------------------------------------------------

# lasso
path.out1 = path.pnsca(N, X = Z, ptype = 1)
boot.out1 = boot.pnsca(N, X = Z, lambda.seq = path.out1$lambda.seq, ptype = 1, trace = TRUE)
out1 = reg.nsca(N, X = Z, penalties = c(boot.out1$lambda.1se, 0, 0), trace = FALSE)

rownames(out1$B) = colnames(Z)
rownames(out1$V) = colnames(N)

round(out1$B, digits = 3)
round(out1$V, digits = 4)

out1$pen.loss
# 0.01001642
out1$penalty
# 0.00209605
out1$SS
# 0.007920372

myplot1 = egg::ggarrange(plots = path.out1$figures, ncol = 2)
ggsave("~/surfdrive/LogitMDA/rcm-association/figures_issp/issp_path.pdf", plot = myplot1, width = 11.7, height = 8.3, units = "in", limitsize = FALSE)

myplot2 = boot.out1$plot
ggsave("~/surfdrive/LogitMDA/rcm-association/figures_issp/issp_boot.pdf", plot = myplot2, width = 11.7, height = 8.3, units = "in", limitsize = FALSE)

#----------------------------------------------------------------------------
# 3 and 4 dimensional analyses
path.out3 = path.pnsca(N, X = Z, S = 3, ptype = 1)
boot.out3 = boot.pnsca(N, X = Z, S = 3, lambda.seq = path.out3$lambda.seq, ptype = 1, trace = TRUE)
out3 = reg.nsca(N, X = Z, S = 3, penalties = c(boot.out3$lambda.1se, 0, 0), trace = FALSE)
rownames(out3$B) = colnames(Z); rownames(out3$V) = colnames(N)
round(out3$B, digits = 3)
out3$pen.loss
# 0.009997022
out3$penalty
# 0.01001642
out3$SS
# 0.007875098

# 4D
path.out4 = path.pnsca(N, X = Z, S = 4, ptype = 1)
boot.out4 = boot.pnsca(N, X = Z, S = 4, lambda.seq = path.out4$lambda.seq, ptype = 1, trace = TRUE)
out4 = reg.nsca(N, X = Z, S = 4, penalties = c(boot.out4$lambda.1se, 0, 0), trace = FALSE)
rownames(out4$B) = colnames(Z); rownames(out4$V) = colnames(N)
round(out4$B, digits = 3)
out4$pen.loss
# 0.009997021
out4$penalty
# 0.002121875
out4$SS
# 0.007875146

#----------------------------------------------------------------------------
# Figure of final solution Constrained NSCA with Lasso penalty
#----------------------------------------------------------------------------

Xo = X
VV = as.data.frame(out1$V)
colnames(VV) = c("dim1", "dim2")
UU = as.data.frame(out1$U)
colnames(UU) = c("dim1", "dim2")
BB = as.data.frame(round(out1$B[-1, ], digits = 3))
colnames(BB) = c("dim1", "dim2")

# lambda scaling - gower e.a. page 24
scaling = (2* nrow(UU)/nrow(VV) / sum(UU^2))^(1/4)

# scaling1 = (nrow(UU)/nrow(VV) * sum(UU[,1]^2))^(1/4)
# scaling2 = (nrow(UU)/nrow(VV) * sum(UU[,2]^2))^(1/4)
# scaling = diag(c(scaling1, scaling2))

UU = UU * scaling
BB = BB * scaling
VV = VV / scaling

# B = BB 
# beta <- B[,2]/B[,1]
# idx.b = which(!is.nan(beta))
# beta = beta[idx.b]
# xnames = rownames(B)[idx.b]
# B = B[idx.b, ]
# Xo = Xo[, idx.b]
# bfin = which(is.finite(beta))
# 
# ynames = colnames(N)
# xynames = c(xnames, ynames)
# 
# beta2 = c(beta, VV[, 2]/VV[, 1])

source("~/surfdrive/LogitMDA/rcm-association/plot2biplot2.R")
myplt = plot2biplot2(UU, VV, BB, Xo)
ggsave("~/surfdrive/LogitMDA/rcm-association/figures_issp/issp_solution3.pdf", plot = myplt, width = 11.7, height = 6.0, units = "in", limitsize = FALSE)



#----------------------------------------------------------------------------
# Predictions based on final 2D model (out1)
#----------------------------------------------------------------------------

PIhat = Z %*% out1$B %*% t(out1$V) + matrix(1, 34, 5) %*% Dc
Nhat = diag(rowSums(N)) %*% PIhat
sum(((Nhat - N)^2/Nhat))

ZX = scale(X)
mx = attr(ZX,"scaled:center")
sdx = attr(ZX,"scaled:scale")

# all predictors at mean, LFPf from 0 to 100 %
Xnew.lfpf = outer(rep(1, 11), colMeans(X)); Xnew.lfpf[, 2] = seq(0, 100, by = 10)
Xnew.lfpf = scale(Xnew.lfpf, center = mx, scale = sdx)
Xnew.lfpf # check - correct
Xnew.lfpf = cbind(1,Xnew.lfpf)
Dc = diag(colSums((N/sum(N))))

# predicted profiles
PIhat = Xnew.lfpf %*% out1$B %*% t(out1$V) + matrix(1, nrow(Xnew.lfpf), 5) %*% Dc
colnames(PIhat) = colnames(N)
round(PIhat, 3)

piplot = data.frame(probs = matrix(PIhat, ncol = 1), 
                    resp = rep(colnames(PIhat), each = 11),
                    x = rep(seq(0,100, by = 10), 5))

ggplot(data = piplot, aes(x = x, y = probs, group=resp)) + geom_line(aes(col= resp)) + geom_point(aes(col= resp)) +
  xlab("Labor Force Participation of Females") + 
  ylab("Predicted Profile")


# [,1]   [,2]   [,3]   [,4]   [,5]
# [1,] 0.560  0.210  0.105  0.015  0.110
# [2,] 0.616  0.180  0.092  0.013  0.099
# [3,] 0.672  0.151  0.078  0.011  0.088
# [4,] 0.728  0.121  0.065  0.010  0.076
# [5,] 0.784  0.092  0.052  0.008  0.065
# [6,] 0.840  0.062  0.038  0.006  0.054
# [7,] 0.896  0.033  0.025  0.004  0.043
# [8,] 0.952  0.003  0.012  0.002  0.031
# [9,] 1.008 -0.026 -0.001  0.000  0.020
# [10,] 1.064 -0.056 -0.015 -0.002  0.009
# [11,] 1.120 -0.085 -0.028 -0.004 -0.003

#-----------------------------------------------------
# blijven niet binnen 0/1 range!
# per rij tellen ze wel op tot 1
#-----------------------------------------------------

# predictions for India with varying levels of LaborForceParticipation of Females (observed value = 21.8%)
Xnew16.lfpf = outer(rep(1, 11), X[16, ]); Xnew16.lfpf[, 2] = seq(0, 100, by = 10)
Xnew16.lfpf = scale(Xnew16.lfpf, center = mx, scale = sdx)
Xnew16.lfpf = cbind(1, Xnew16.lfpf) 

# predicted profiles
PIhat = Xnew16.lfpf %*% out1$B %*% t(out1$V) + matrix(1, nrow(Xnew16.lfpf), 5) %*% Dc
colnames(PIhat) = colnames(N)
round(PIhat, 3)

piplot = data.frame(probs = matrix(PIhat, ncol = 1), 
                    resp = rep(colnames(PIhat), each = 11),
                    x = rep(seq(0,100, by = 10), 5))

piplot$probs = pmin(piplot$probs, 1)
piplot$probs = pmax(piplot$probs, 0)

plt1 = ggplot(data = piplot, aes(x = x, y = probs, group= resp)) + 
  geom_line(aes(col= resp)) + geom_point(aes(col= resp)) +
  xlab("Labor Force Participation of Females") + 
  ylab("Predicted Profile") + 
  scale_x_continuous(breaks=seq(0, 100, 10)) + 
  scale_y_continuous(breaks=seq(-0.1, 1.1, .1)) + 
  theme_bw()

# ##################################################################
# # variable axes with markers
# ##################################################################
# # for solid line
# MCx1 <- data.frame(labs=character(),
#                    varx = integer(),
#                    dim1 = double(),
#                    dim2 = double(), stringsAsFactors=FALSE)
# # for markers
# MCx2 <- data.frame(labs=character(),
#                    varx = integer(),
#                    dim1 = double(),
#                    dim2 = double(), stringsAsFactors=FALSE)
# 
# ll = 0
# lll = 0
# for(p in 1:ncol(Xo)){
#   b = matrix(as.numeric(B[p , ]), 2, 1)
#   # solid line
#   minx = min(Xo[, p])
#   maxx = max(Xo[, p])
#   m.x1 = c(minx,maxx)
#   markers1 = matrix((m.x1 - mx[p])/sdx[p], 2, 1)
#   markerscoord1 = outer(markers1, b) # markers1 %*% t(b %*% solve(t(b) %*% b))
#   MCx1[(ll + 1): (ll + 2), 1] = paste0(c("min", "max"), p)
#   MCx1[(ll + 1): (ll + 2), 2] = p
#   MCx1[(ll + 1): (ll + 2), 3:4] = markerscoord1
#   ll = ll + 2
#   # markers
#   m.x2 = pretty(Xo[, p])
#   m.x2 = m.x2[which(m.x2 > minx & m.x2 < maxx)]
#   l.m = length(m.x2)
#   markers2 = matrix((m.x2 - mx[p])/sdx[p], l.m, 1)
#   markerscoord2 = outer(markers2, b) # markers2 %*% t(b %*% solve(t(b) %*% b))
#   MCx2[(lll + 1): (lll + l.m), 1] = paste(m.x2)
#   MCx2[(lll + 1): (lll + l.m), 2] = p
#   MCx2[(lll + 1): (lll + l.m), 3:4] = markerscoord2
#   lll = lll + l.m
# } # loop p
# 
# ##################################################################
# # variable axes with markers
# ##################################################################
# # for solid line
# MCy1 <- data.frame(labs=character(),
#                    varx = integer(),
#                    dim1 = double(),
#                    dim2 = double(), stringsAsFactors=FALSE)
# # for markers
# MCy2 <- data.frame(labs=character(),
#                    varx = integer(),
#                    dim1 = double(),
#                    dim2 = double(), stringsAsFactors=FALSE)
# 
# PI = out1$PI
# ll = 0
# lll = 0
# for(p in 1:ncol(N)){
#   v = matrix(as.numeric(VV[p , ]), 2, 1)
#   # solid line
#   minx = min(PI[, p])
#   maxx = max(PI[, p])
#   m.x1 = c(minx,maxx)
#   markers1 = matrix(m.x1, 2, 1)
#   markerscoord1 = markers1 %*% t(v %*% solve(t(v) %*% v))
#   MCy1[(ll + 1): (ll + 2), 1] = paste0(c("min", "max"), p)
#   MCy1[(ll + 1): (ll + 2), 2] = p
#   MCy1[(ll + 1): (ll + 2), 3:4] = markerscoord1
#   ll = ll + 2
#   # markers
#   m.x2 = pretty(PI[, p])
#   m.x2 = m.x2[which(m.x2 > minx & m.x2 < maxx)]
#   l.m = length(m.x2)
#   markers2 = matrix(m.x2, l.m, 1)
#   markerscoord2 = markers2 %*% t(v %*% solve(t(v) %*% v)) # outer(markers2, v)
#   MCy2[(lll + 1): (lll + l.m), 1] = paste(m.x2)
#   MCy2[(lll + 1): (lll + l.m), 2] = p
#   MCy2[(lll + 1): (lll + l.m), 3:4] = markerscoord2
#   lll = lll + l.m
# } # loop p
# 
# 
# xcol = "lightskyblue"
# ycol = "red"
# 
# plt = ggplot() + 
# #  geom_text(data = VV, aes(x = dim1, y = dim2, label = rownames(VV)), colour = "red")  + 
#   geom_text(data = UU, aes(x = dim1, y = dim2, label = rownames(UU)), colour = "blue")
# 
# margins <- c("l" = ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range[1] - .1,
#              "r" = ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range[2] + .1,
#              "b" = ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range[1] - .1,
#              "t" = ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range[2] + .1)
# 
# if(any(!is.finite(beta))) plt = plt + geom_vline(xintercept = 0, colour = xcol, linetype = 3)
# 
# plt = plt + geom_abline(intercept = 0, slope = beta[bfin], colour = xcol, linetype = 3) +
#   geom_line(data = MCx1, aes(x = dim1, y = dim2, group = varx), col = xcol, linewidth = 1) +
#   geom_point(data = MCx2, aes(x = dim1, y = dim2), col = xcol) +
#   geom_text(data = MCx2, aes(x = dim1, y = dim2, label = labs), nudge_y = -0.01, size = 1.5)
# 
# plt = plt + geom_abline(intercept = 0, slope = (VV[, 2]/VV[, 1]), colour = ycol, linetype = 3)+
#   geom_line(data = MCy1, aes(x = dim1, y = dim2, group = varx), col = ycol, linewidth = 1) +
#   geom_point(data = MCy2, aes(x = dim1, y = dim2), col = ycol) +
#   geom_text(data = MCy2, aes(x = dim1, y = dim2, label = labs), nudge_y = -0.01, size = 1.5)
#   
# 
# ######################################################
# # add x-variable labels
# ######################################################
# 
# lab <- data.frame("xname" = xynames, 
#                   "b" = beta2, 
#                   "Yleft" = beta2*margins["l"],
#                   "Yright" = beta2*margins["r"])
# 
# BBB = rbind((B + (B==0) * 1e6), VV) # trick to display name of GDP
# BV = rbind(B, VV)
# orientation = sign(BBB[,1]) #sign of dim1 defines direction l-r
# lab$side =  c("left","right")[ as.numeric(BV[,1] > 0)+1] 
# lab$side[lab$Yleft < margins["b"] & orientation<0 ] = "bottom"
# lab$side[lab$Yleft > margins["t"] & orientation<0 ] = "top"
# lab$side[lab$Yright < margins["b"]& orientation>0] = "bottom"
# lab$side[lab$Yright > margins["t"]& orientation>0] = "top"
# 
# lab$X <- lab$Y <- NA
# lab$X[lab$side == "bottom"] <- (margins["b"]/beta2[lab$side == "bottom"])
# lab$X[lab$side == "top"] <- (margins["t"]/beta2[lab$side == "top"])
# lab$Y[lab$side == "left"] <- margins["l"]*beta2[lab$side == "left"]
# lab$Y[lab$side == "right"] <-margins["r"]*beta2[lab$side == "right"]
# 
# lab <- split(lab, lab$side)
# 
# plt = plt + 
#   scale_x_continuous(breaks = lab$bottom$X, labels = lab$bottom$xname, sec.axis = sec_axis(trans ~ ., breaks = lab$top$X, labels = lab$top$xname)) +
#   scale_y_continuous(breaks = lab$left$Y, labels = lab$left$xname, sec.axis = sec_axis(trans ~ ., breaks = lab$right$Y, labels = lab$right$xname))
# 
# myplot3 = plt + coord_fixed(xlim = margins[c("l","r")], ylim = margins[c("b","t")], expand = F) +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())
# 
# myplot3
# ggsave("~/surfdrive/LogitMDA/rcm-association/figures_issp/issp_solution.pdf", plot = myplot3, width = 11.7, height = 8.3, units = "in", limitsize = FALSE)
