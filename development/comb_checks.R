devtools::load_all('.', quiet=TRUE)
data(caldana)

cat('=== comb() basic test ===\n')
mod1 <- tryCatch(
  asca(compounds ~ comb(light + time), data=caldana),
  error = function(e) { cat('FAIL:', conditionMessage(e), '\n'); NULL }
)
if (!is.null(mod1)) cat('comb(light+time): OK, effects=', paste(names(mod1$ssq), collapse=', '), '\n')

cat('\n=== comb() + main effect ===\n')
mod2 <- tryCatch(
  asca(compounds ~ time + comb(light + time:light), data=caldana),
  error = function(e) { cat('FAIL:', conditionMessage(e), '\n'); NULL }
)
if (!is.null(mod2)) cat('time+comb(light+time:light): OK, effects=', paste(names(mod2$ssq), collapse=', '), '\n')

cat('\n=== comb() + permutation ===\n')
if (!is.null(mod1)) {
  p <- tryCatch(
    permutation(mod1, permute=50),
    error = function(e) { cat('FAIL perm:', conditionMessage(e), '\n'); NULL }
  )
  if (!is.null(p)) cat('permutation on comb(): OK, pvalues=', round(p$permute$pvalues, 3), '\n')
}

cat('\n=== comb() respect_SStype ===\n')
mod3 <- tryCatch(
  asca(compounds ~ comb(light + time), data=caldana, respect_SStype=TRUE),
  error = function(e) { cat('FAIL:', conditionMessage(e), '\n'); NULL }
)
if (!is.null(mod3)) cat('comb() respect_SStype=TRUE: OK\n')

cat('\n=== comb() SStype I/II/III ===\n')
for (ss in 1:3) {
  m <- tryCatch(
    asca(compounds ~ comb(light + time), data=caldana, SStype=ss),
    error = function(e) { cat('FAIL SStype', ss, ':', conditionMessage(e), '\n'); NULL }
  )
  if (!is.null(m)) cat('SStype', ss, ': OK\n')
}

cat('\nAll comb() checks done.\n')
