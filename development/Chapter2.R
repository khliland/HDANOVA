library(latex2exp)
X <- matrix(c(3,1,5,1,3,2), 3,2)

# Filled circles with aquamarine-like color
x11()
par(mar=c(4,4,1,1), mgp=c(2,0.7,0))
plot(X[,1],X[,2], cex = 3, col = rgb(0,0.5,0.5,0.5), pch = 19,
     xlab = TeX(r'($V_1$)'), ylab = TeX(r'($V_2$)'),
     xlim = c(0.8, 5.2), ylim = c(0.8, 3.2))
for(i in 1:3){
  text(X[i,1]-0.05,X[i,2]+0.1, labels = TeX(sprintf(r'($O_%d$)',i)), pos = 2)
}

# Three dimensional plot with filled circles, where the rows of X are the three coordinates values for each of the two observations
# Filled circles with aquamarine-like color
s3d <- scatterplot3d(t(X), xlim = c(0,5), ylim = c(0,5), zlim = c(0,5), pch = 19, color = rgb(0,0.5,0.5,0.5), cex.symbols = 3,
              xlab = TeX(r'($O_1$)'), ylab = TeX(r'($O_2$)'), zlab = TeX(r'($O_3$)'), box = FALSE, grid = FALSE)

a <- s3d$xyz.convert(3,1,5); b <- s3d$xyz.convert(1,3,2)
text(a$x-0.05, a$y+0.1, labels = TeX(r'($V_1$)'), pos = 2)
text(b$x-0.05, b$y+0.1, labels = TeX(r'($V_2$)'), pos = 2)

# Add support lines
a <- s3d$xyz.convert(3,1,5); b <- s3d$xyz.convert(3,1,0)
lines(c(a$x, b$x), c(a$y, b$y), col = "gray", lwd = 2)
a <- s3d$xyz.convert(3,1,0); b <- s3d$xyz.convert(3,0,0)
lines(c(a$x, b$x), c(a$y, b$y), col = "gray", lwd = 2)
a <- s3d$xyz.convert(3,1,0); b <- s3d$xyz.convert(5,1,0)
lines(c(a$x, b$x), c(a$y, b$y), col = "gray", lwd = 2)

a <- s3d$xyz.convert(1,3,2); b <- s3d$xyz.convert(1,3,0)
lines(c(a$x, b$x), c(a$y, b$y), col = "gray", lwd = 2)
a <- s3d$xyz.convert(1,3,0); b <- s3d$xyz.convert(1,0,0)
lines(c(a$x, b$x), c(a$y, b$y), col = "gray", lwd = 2)
a <- s3d$xyz.convert(1,3,0); b <- s3d$xyz.convert(5,3,0)
lines(c(a$x, b$x), c(a$y, b$y), col = "gray", lwd = 2)

