food = matrix(
  c(0, 0, 2, 0, 0, 15,
    0, 1, 0, 22, 11, 9,
    2, 17, 0, 0, 5, 9,
    0, 4, 0, 17, 12, 6,
    4, 3, 0, 0, 1, 7,
    0, 0, 0, 3, 14, 0,
    0, 10, 0, 0, 11, 0,
    3, 16, 5, 10, 16, 13,
    8, 9, 0, 2, 0, 0,
    1, 2, 6, 5, 2, 12,
    10, 4, 1, 5, 2, 9,
    14, 3, 0, 3, 0, 25, 
    4, 1, 3, 2, 11, 6, 
    8, 0, 0, 0, 12, 1, 
    1, 1, 0, 0, 5, 12,
    0, 0, 8, 2, 12, 9,
    4, 0, 2, 1, 4, 2, 
    1, 6, 1, 0, 11, 2, 
    2, 2, 5, 2, 5, 5, 
    2, 8, 0, 6, 11, 22,
    10, 0, 0, 5, 0, 0, 
    3, 4, 0, 1, 0, 4,
    0, 7, 0, 0, 6, 0,
    1, 7, 4, 2, 1, 1, 
    5, 5, 1, 2, 1, 1,
    0, 12, 1, 19, 5, 11,
    2, 1, 3, 1, 3, 1, 
    8, 3, 4, 0, 17, 8), 6, 28)
rownames(food) = c("Belgium", "France", "Italy", "Norway", "Poland", "Spain")
colnames(food) = c("Ancient", "Christmas", "Cooking", "Country", "Culture", "Diner", "Dish", "Family",
                   "Feast", "Good", "Grandmother", "Habit", "Healthy", "Holidays", "Home", "Home-made",
                   "Kitchen", "Meal", "Natural", "Old", "Old-fashioned", "Quality", "Recipe", 
                   "Regional", "Restaurant", "Rural", "Simple", "Tasty")
colnames(food) = as.character(seq(1:28))
food
xtable::xtable(food)

X = diag(6)
colnames(X) = rownames(food)

library(CAvariants)
out = CAvariants(t(food), catype = "NSCA")
summary(out)

# Total inertia  0.033 
# 
# Inertias, percent inertias and cumulative percent inertias of the row and column space
# 
# inertia inertiapc cuminertiapc
# value1   0.011    32.401       32.401
# value2   0.008    23.828       56.229
# value3   0.006    18.735       74.964
# value4   0.005    15.604       90.568
# value5   0.003     9.432      100.000
# 
# Predictability Index for Variants of Non symmetrical Correspondence Analysis:
#   
#   Tau Index predicting from column 
# 
# [1] 0.035
# 
# C-statistic 698.342 and p-value 0

source("~/surfdrive/LogitMDA/rcm-association/nsca.R")

path.out = path.pnsca(food, X = X, ptype = 3)
boot.out = boot.pnsca(food, X = X, lambda.seq = path.out$lambda.seq, ptype = 3, trace = TRUE)
out = reg.nsca(food, X = X, penalties = c(0, 0, boot.out$lambda.min), trace = FALSE)

rownames(out$B) = rownames(out$U) = rownames(food)
rownames(out$V) = c("Ancient", "Christmas", "Cooking", "Country", "Culture", "Diner", "Dish", "Family",
                    "Feast", "Good", "Grandmother", "Habit", "Healthy", "Holidays", "Home", "Home-made",
                    "Kitchen", "Meal", "Natural", "Old", "Old-fashioned", "Quality", "Recipe", 
                    "Regional", "Restaurant", "Rural", "Simple", "Tasty")
round(out$B, digits = 3)
round(out$V, digits = 3)

# explained variance
out$pen.loss
out$SS
(0.033 - out$SS)/0.033

# penalty value
out$pen.loss - out$SS
boot.out$lambda.min * sum(sqrt(rowSums(out$B^2)))
boot.out$lambda.1se * sum(sqrt(rowSums(out$B^2)))

# # save figures
# p1 = egg::ggarrange(plots = path.out$figures, ncol = 2)
# ggsave("~/surfdrive/LogitMDA/rcm-association/figures_food/food_path.pdf", plot = p1, width = 11.7, height = 8.3, units = "in", limitsize = FALSE)
# ggsave("~/surfdrive/LogitMDA/rcm-association/figures_food/food_boot.pdf", plot = boot.out$plot, width = 11.7, height = 8.3, units = "in", limitsize = FALSE)
# 

#----------------------------------------------------------------------------
# Figure of final solution NSCA with Group Lasso penalty
#----------------------------------------------------------------------------

# lambda scaling - gower e.a. page 24
scaling = (2* nrow(out$U)/nrow(out$V) / sum(out$U^2))^(1/4)

VV = as.data.frame(out$V / scaling)
colnames(VV) = c("dim1", "dim2")
UU = as.data.frame(out$U * scaling)
colnames(UU) = c("dim1", "dim2")

ynames = colnames(food)
beta2 = VV[ , 2]/VV[, 1]
  
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

PI = out$PI
ll = 0
lll = 0
for(p in 1:ncol(food)){
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

plt = ggplot() + geom_blank(data = VV, aes(x = dim1, y = dim2))

margins <- c("l" = ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range[1] - .1,
             "r" = ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range[2] + .1,
             "b" = ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range[1] - .1,
             "t" = ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range[2] + .1)

plt = plt + geom_abline(intercept = 0, slope = (VV[, 2]/VV[, 1]), colour = ycol, linetype = 3)+
  geom_line(data = MCy1, aes(x = dim1, y = dim2, group = varx), col = ycol) +
  geom_point(data = MCy2, aes(x = dim1, y = dim2), col = ycol) +
  geom_text(data = MCy2, aes(x = dim1, y = dim2, label = labs), nudge_y = -0.01, size = 1.5)

UUU = UU[4:5, ] #lambda_1se
UUU = UU[c(1,4,5,6), ] # lambda_min

plt = plt + 
  geom_point(data = UU, aes(x = dim1, y = dim2), size = 4, colour = "blue") + 
  geom_text(data = UUU, aes(x = dim1, y = dim2, label = rownames(UUU)), colour = "blue", nudge_y = -0.02, size = 3)


######################################################
# add x-variable labels
######################################################


lab <- data.frame("xname" = ynames, 
                  "b" = beta2, 
                  "Yleft" = beta2*margins["l"],
                  "Yright" = beta2*margins["r"])

orientation = sign(VV[,1]) #sign of dim1 defines direction l-r
lab$side =  c("left","right")[ as.numeric(VV[,1] > 0)+1] 
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
ggsave("~/surfdrive/LogitMDA/rcm-association/figures_food/food_solution.pdf", plot = myplot3, width = 11.7, height = 8.3, units = "in", limitsize = FALSE)
