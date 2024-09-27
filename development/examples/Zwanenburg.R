x = c(seq(0.25,1,by=0.25), seq(1.5,6,by=0.5))
alpha1 = c(0,0,0,0,0,0.005,0.027,0.085, 0.172, 0.228, 0.26, 0.322, 0.35,0.372)
alpha2 = c(0,0,0,0,0,0.0125,0.1,0.19, 0.25, 0.345, 0.36, 0.4, 0.411,0.445)
beta1 = c(0,0,0.005,0.04, 0.21, 0.34, 0.375, 0.422, 0.4375, 0.44,0.46,0.465,0.475, 0.47)
beta2 = c(0,0,0.022,0.128, 0.3125, 0.4, 0.43, 0.438, 0.47, 0.47,0.48,0.48,0.5, 0.488)

pdf("Figure_4.16.pdf", 18,8, pointsize = 18)
par(mfrow=c(1,2), mar=c(3,4,3,1), mgp=c(3,0.5,0), xpd=NA, tcl=0.5)
# Set the xlim so that the plot starts exactly at 0 and ends exactly at 6, including the 5%
scale = 0.037
matplot(x, cbind(alpha1, alpha2), type="l", lty=1, lwd=2, col=1, ylab="p-value",
        cex=0.5, main=expression(paste("Factor ", alpha)),
        xlab="",
        # Compensate xlim and ylim such that usr will be c(0,6,0,0.45):
        xlim = c(0+6*scale,6-6*scale), ylim=c(0+0.45*scale, 0.45-0.45*scale),
        axes = FALSE)
axis(1, at=seq(0,6,by=1))
axis(2, at=seq(0,0.45,by=0.05), las=2)
box()
points(x, alpha1, col=1, pch=8, cex=0.9)
points(x, alpha2, col=1, pch=1, cex=1.1)
mtext(expression(paste("Noise level ", delta)), 1, 2)

matplot(x, cbind(beta1, beta2), type="l", lty=1, lwd=2, col=1, ylab="p-value",
        cex=0.5, main=expression(paste("Factor ", beta)),
        xlab="",
        xlim = c(0+6*scale,6-6*scale), ylim=c(0+0.5*scale, 0.5-0.5*scale),axes = FALSE)
axis(1, at=seq(0,6,by=1))
axis(2, at=seq(0,0.5,by=0.05), las=2)
box()
points(x, beta1, col=1, pch=8, cex=0.9)
points(x, beta2, col=1, pch=1, cex=1.1)
mtext(expression(paste("Noise level ", delta)), 1, 2)
dev.off()
