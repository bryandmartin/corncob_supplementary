## @knitr metrics
library(corncob)
n <- new_metric("n", "n", 
                metric = function(model, out) {
                  length(out$mod$W)
})

cons <- new_metric("cons", "constant", 
                metric = function(model, out) {
                  out$cons
                })

wald <- new_metric("wald", "Wald Test",
                        metric = function(model, out) {
                          out$wald
})

pbwald <- new_metric("pbwald", "Parametric Bootstrap Wald Test",
                     metric = function(model, out) {
                       out$pbwald
                     })

lrt <- new_metric("lrt", "Likelihood Ratio Test",
                  metric = function(model, out) {
                    out$lrt
                  })

pblrt <- new_metric("pblrt", "Parametric Bootstrap Likelihood Ratio Test",
                     metric = function(model, out) {
                       out$pblrt
                     })

