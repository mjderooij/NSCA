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

#----------------------------------------------------------------------------
# Constrained NSCA with Lasso penalty
#----------------------------------------------------------------------------

# lasso
path.out1 = path.pnsca(N, X = Z, ptype = 1)
boot.out1 = boot.pnsca(N, X = Z, lambda.seq = path.out1$lambda.seq, ptype = 1, trace = TRUE)
out1 = reg.nsca(N, X = Z, penalties = c(boot.out1$lambda.1se, 0, 0), trace = FALSE)

rownames(out1$B) = colnames(Z)
rownames(out1$V) = colnames(N)

round(out1$B, digits = 2)
round(out1$V, digits = 2)


myplot1 = egg::ggarrange(plots = path.out1$figures, ncol = 2)
ggsave("~/surfdrive/LogitMDA/rcm-association/figures_issp/issp_path.pdf", plot = myplot1, width = 11.7, height = 8.3, units = "in", limitsize = FALSE)

myplot2 = boot.out1$plot
ggsave("~/surfdrive/LogitMDA/rcm-association/figures_issp/issp_boot.pdf", plot = myplot2, width = 11.7, height = 8.3, units = "in", limitsize = FALSE)

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

B = BB 
beta <- B[,2]/B[,1]
idx.b = which(!is.nan(beta))
beta = beta[idx.b]
xnames = rownames(B)[idx.b]
B = B[idx.b, ]
Xo = Xo[, idx.b]
bfin = which(is.finite(beta))

ynames = colnames(N)
xynames = c(xnames, ynames)

beta2 = c(beta, VV[, 2]/VV[, 1])

##################################################################
# variable axes with markers
##################################################################
# for solid line
MCx1 <- data.frame(labs=character(),
                   varx = integer(),
                   dim1 = double(),
                   dim2 = double(), stringsAsFactors=FALSE)
# for markers
MCx2 <- data.frame(labs=character(),
                   varx = integer(),
                   dim1 = double(),
                   dim2 = double(), stringsAsFactors=FALSE)

ll = 0
lll = 0
for(p in 1:ncol(Xo)){
  b = matrix(as.numeric(B[p , ]), 2, 1)
  # solid line
  minx = min(Xo[, p])
  maxx = max(Xo[, p])
  m.x1 = c(minx,maxx)
  markers1 = matrix((m.x1 - mx[p])/sdx[p], 2, 1)
  markerscoord1 = outer(markers1, b) # markers1 %*% t(b %*% solve(t(b) %*% b))
  MCx1[(ll + 1): (ll + 2), 1] = paste0(c("min", "max"), p)
  MCx1[(ll + 1): (ll + 2), 2] = p
  MCx1[(ll + 1): (ll + 2), 3:4] = markerscoord1
  ll = ll + 2
  # markers
  m.x2 = pretty(Xo[, p])
  m.x2 = m.x2[which(m.x2 > minx & m.x2 < maxx)]
  l.m = length(m.x2)
  markers2 = matrix((m.x2 - mx[p])/sdx[p], l.m, 1)
  markerscoord2 = outer(markers2, b) # markers2 %*% t(b %*% solve(t(b) %*% b))
  MCx2[(lll + 1): (lll + l.m), 1] = paste(m.x2)
  MCx2[(lll + 1): (lll + l.m), 2] = p
  MCx2[(lll + 1): (lll + l.m), 3:4] = markerscoord2
  lll = lll + l.m
} # loop p

##################################################################
# variable axes with markers
##################################################################
# for solid line
MCy1 <- data.frame(labs=character(),
                   varx = integer(),
                   dim1 = double(),
                   dim2 = double(), stringsAsFactors=FALSE)
# for markers
MCy2 <- data.frame(labs=character(),
                   varx = integer(),
                   dim1 = double(),
                   dim2 = double(), stringsAsFactors=FALSE)

PI = out1$PI
ll = 0
lll = 0
for(p in 1:ncol(N)){
  v = matrix(as.numeric(VV[p , ]), 2, 1)
  # solid line
  minx = min(PI[, p])
  maxx = max(PI[, p])
  m.x1 = c(minx,maxx)
  markers1 = matrix(m.x1, 2, 1)
  markerscoord1 = markers1 %*% t(v %*% solve(t(v) %*% v))
  MCy1[(ll + 1): (ll + 2), 1] = paste0(c("min", "max"), p)
  MCy1[(ll + 1): (ll + 2), 2] = p
  MCy1[(ll + 1): (ll + 2), 3:4] = markerscoord1
  ll = ll + 2
  # markers
  m.x2 = pretty(PI[, p])
  m.x2 = m.x2[which(m.x2 > minx & m.x2 < maxx)]
  l.m = length(m.x2)
  markers2 = matrix(m.x2, l.m, 1)
  markerscoord2 = markers2 %*% t(v %*% solve(t(v) %*% v)) # outer(markers2, v)
  MCy2[(lll + 1): (lll + l.m), 1] = paste(m.x2)
  MCy2[(lll + 1): (lll + l.m), 2] = p
  MCy2[(lll + 1): (lll + l.m), 3:4] = markerscoord2
  lll = lll + l.m
} # loop p


xcol = "lightskyblue"
ycol = "red"

plt = ggplot() + 
#  geom_text(data = VV, aes(x = dim1, y = dim2, label = rownames(VV)), colour = "red")  + 
  geom_text(data = UU, aes(x = dim1, y = dim2, label = rownames(UU)), colour = "blue")

margins <- c("l" = ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range[1] - .1,
             "r" = ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range[2] + .1,
             "b" = ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range[1] - .1,
             "t" = ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range[2] + .1)

if(any(!is.finite(beta))) plt = plt + geom_vline(xintercept = 0, colour = xcol, linetype = 3)

plt = plt + geom_abline(intercept = 0, slope = beta[bfin], colour = xcol, linetype = 3) +
  geom_line(data = MCx1, aes(x = dim1, y = dim2, group = varx), col = xcol, linewidth = 1) +
  geom_point(data = MCx2, aes(x = dim1, y = dim2), col = xcol) +
  geom_text(data = MCx2, aes(x = dim1, y = dim2, label = labs), nudge_y = -0.01, size = 1.5)

plt = plt + geom_abline(intercept = 0, slope = (VV[, 2]/VV[, 1]), colour = ycol, linetype = 3)+
  geom_line(data = MCy1, aes(x = dim1, y = dim2, group = varx), col = ycol, linewidth = 1) +
  geom_point(data = MCy2, aes(x = dim1, y = dim2), col = ycol) +
  geom_text(data = MCy2, aes(x = dim1, y = dim2, label = labs), nudge_y = -0.01, size = 1.5)
  

######################################################
# add x-variable labels
######################################################



lab <- data.frame("xname" = xynames, 
                  "b" = beta2, 
                  "Yleft" = beta2*margins["l"],
                  "Yright" = beta2*margins["r"])

BBB = rbind((B + (B==0) * 1e6), VV) # trick to display name of GDP
BV = rbind(B, VV)
orientation = sign(BBB[,1]) #sign of dim1 defines direction l-r
lab$side =  c("left","right")[ as.numeric(BV[,1] > 0)+1] 
lab$side[lab$Yleft < margins["b"] & orientation<0 ] = "bottom"
lab$side[lab$Yleft > margins["t"] & orientation<0 ] = "top"
lab$side[lab$Yright < margins["b"]& orientation>0] = "bottom"
lab$side[lab$Yright > margins["t"]& orientation>0] = "top"

lab$X <- lab$Y <- NA
lab$X[lab$side == "bottom"] <- (margins["b"]/beta2[lab$side == "bottom"])
lab$X[lab$side == "top"] <- (margins["t"]/beta2[lab$side == "top"])
lab$Y[lab$side == "left"] <- margins["l"]*beta2[lab$side == "left"]
lab$Y[lab$side == "right"] <-margins["r"]*beta2[lab$side == "right"]

lab <- split(lab, lab$side)

plt = plt + 
  scale_x_continuous(breaks = lab$bottom$X, labels = lab$bottom$xname, sec.axis = sec_axis(trans ~ ., breaks = lab$top$X, labels = lab$top$xname)) +
  scale_y_continuous(breaks = lab$left$Y, labels = lab$left$xname, sec.axis = sec_axis(trans ~ ., breaks = lab$right$Y, labels = lab$right$xname))

myplot3 = plt + coord_fixed(xlim = margins[c("l","r")], ylim = margins[c("b","t")], expand = F) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

myplot3
ggsave("~/surfdrive/LogitMDA/rcm-association/figures_issp/issp_solution.pdf", plot = myplot3, width = 11.7, height = 8.3, units = "in", limitsize = FALSE)
