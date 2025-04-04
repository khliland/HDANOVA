## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=6, 
  fig.height=4
)
# Legge denne i YAML på toppen for å skrive ut til tex
#output: 
#  pdf_document: 
#    keep_tex: true
# Original:
#  rmarkdown::html_vignette:
#    toc: true

## ----setup--------------------------------------------------------------------
# Start the mixlm R package
library(mixlm, warn.conflicts = FALSE)

## -----------------------------------------------------------------------------
set.seed(123)
dat <- data.frame(
  feed  = factor(rep(rep(c("low","high"), each=6), 4)),
  breed = factor(rep(c("NRF","Hereford","Angus"), 16)),
  bull  = factor(rep(LETTERS[1:4], each = 12)),
  daughter = factor(c(rep(letters[1:4], 3), rep(letters[5:8], 3), rep(letters[9:12], 3), rep(letters[13:16], 3))),
  age   = round(rnorm(48, mean = 36, sd = 5))
)
dat$yield <- 150*with(dat, 10 + 3 * as.numeric(feed) + as.numeric(breed) + 
                        2 * as.numeric(bull) + 1 * as.numeric(sample(dat$daughter, 48)) + 
                        0.5 * age + rnorm(48, sd = 2))
head(dat)

## -----------------------------------------------------------------------------
mod <- lm(yield ~ feed, data = dat)
print(anova(mod))

## -----------------------------------------------------------------------------
summary(mod)

## ----fig.width=6, fig.height=5------------------------------------------------
old.par <- par(mfrow=c(2,1), mar=c(4,4,2,0.5))
plot(mod, which = 1:2, ask=FALSE)
par(old.par)

## -----------------------------------------------------------------------------
mod <- lm(yield ~ feed*breed, data = dat)
print(anova(mod))

## -----------------------------------------------------------------------------
mod <- lm(yield ~ feed+breed, data = dat)
print(anova(mod))

## -----------------------------------------------------------------------------
# Type I - Sequential testing, including one and one effect
print(anova(mod))

# Type II - Testing each term after all others, 
# except ignoring the term's higher-order relatives
print(Anova(mod, type="II"))

# Type III - Testing each term after all others,
# including the term's higher-order relatives
print(Anova(mod, type="III"))

## -----------------------------------------------------------------------------
# Sum-coding, i.e., the sum of all levels is zero and all effects
# are orthogonal in the balanced case.
mod <- lm(yield ~ feed*breed, data = dat, contrasts="contr.sum")
print(Anova(mod, type="III"))

# Weighted coding, i.e., the sum of all levels is zero and the effects
# are weighted by the number of levels, effect-wise.
mod <- lm(yield ~ feed*breed, data = dat, contrasts="contr.weighted")
print(Anova(mod, type="III"))

## -----------------------------------------------------------------------------
options(contrasts = c("contr.sum", "contr.poly"))

## -----------------------------------------------------------------------------
mod <- lm(yield ~ feed*breed + age, data = dat)
print(Anova(mod, type="II"))

## -----------------------------------------------------------------------------
mod <- lm(yield ~ bull + daughter%in%bull, data = dat)
print(Anova(mod, type="II"))

## -----------------------------------------------------------------------------
mod <- lm(yield ~ feed*r(bull), data = dat)
print(Anova(mod, type="II"))

## -----------------------------------------------------------------------------
mod <- lm(yield ~ feed*r(bull), data = dat, unrestricted = FALSE)
print(Anova(mod, type="II"))

## -----------------------------------------------------------------------------
set.seed(123)
long <- dat[c(1:4,9:12), c("feed", "daughter", "yield")]
long <- rbind(long, long, long)
long$daughter <- factor(long$daughter) # Remove redundant daughter
long$time  <- factor(rep(1:3, each=8))
long$yield <- long$yield + rnorm(24, sd = 100) + rep(c(-200,0,200), each=8)
plot(yield~daughter, data=long)

## -----------------------------------------------------------------------------
mod <- lm(yield ~ r(daughter) + feed*r(time), data = long, unrestricted=FALSE)
print(Anova(mod, type="II"))

## -----------------------------------------------------------------------------
mod <- lm(yield ~ feed*r(bull), data = dat, REML = TRUE)
print(Anova(mod, type="III"))

## -----------------------------------------------------------------------------
print(mod)

## -----------------------------------------------------------------------------
dat$mastitis <- as.numeric(dat$breed) + as.numeric(dat$feed) + rnorm(48, sd = 1)

## -----------------------------------------------------------------------------
mod <- lm(cbind(yield,mastitis) ~ feed*breed, data = dat)
print(Anova(mod, type="II"))

## -----------------------------------------------------------------------------
print(Anova(mod, type="II", test="Wilks"))
print(Anova(mod, type="II", test="Hotelling-Lawley"))
print(Anova(mod, type="II", test="Roy"))

## -----------------------------------------------------------------------------
summary(mod)

