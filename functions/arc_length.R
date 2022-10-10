library(pracma)
library(TestGardener)
library(mirt)


arc_length <- function(mod, from, to, surp = T) {
  if (class(mod)=="SingleGroupClass") { # check if mirt model
    if (surp) arclength(mirt_surps, from, to, mirt_mod = mod)
    else arclength(mirt_probs, from, to, mirt_mod = mod)
  }
  else { # assume testGardener
    if (surp) arclength(tg_surps, from, to, analyze_res = mod)
    else arclength(tg_probs, from, to, analyze_res = mod)
  }
}

mirt_tg_arc_length <- function(mirt_mod, tg_mod, mirt_from, mirt_to, 
                               mirt_tg_lm, surp = T) {
  if (surp) arclength(mirt_tg_surps, from, to, 
                      mirt_mod = mod, tg_mod = tg_mod,
                      mirt_thetas = mirt_thetas, tg_thetas = tg_thetas)
  else arclength(mirt_tg_probs, from, to, mirt_mod = mod)
}

# takes thetas from mirt and equivalent from thetas from TG as well as the models
# returns probabilities for each response category
mirt_tg_probs <- function(mirt_theta, 
                          mirt_mod, tg_mod, 
                          mirt_thetas, tg_thetas) {
  mirt_probs <- c(probtrace(mirt_mod, mirt_theta))
  
  tg_theta = approx(mirt_thetas, tg_thetas, xout = mirt_theta)
  
  iter <- length(analyze_res$parList)
  no_items <- length(analyze_res$parList[[iter]]$WfdList)
  tg_probs <- vector("numeric", 0)
  for (item in 1:no_items) {
    WListi <- analyze_res$parList[[iter]]$WfdList[[item]]
    Wfdi <- WListi$Wfd
    Mi <- WListi$M
    logMi <- log(Mi)
    surpisal <- eval.surp(theta, Wfdi)
    tg_probs <- c(tg_probs, exp(-surpisal * logMi))
  }
  c(mirt_probs, tg_probs)
}

tg_probs <- function(theta, analyze_res) {
  iter <- length(analyze_res$parList)
  no_items <- length(analyze_res$parList[[iter]]$WfdList)
  probs <- vector("numeric", 0)
  for (item in 1:no_items) {
    WListi <- analyze_res$parList[[iter]]$WfdList[[item]]
    Wfdi <- WListi$Wfd
    Mi <- WListi$M
    logMi <- log(Mi)
    surpisal <- eval.surp(theta, Wfdi)
    probs <- c(probs, exp(-surpisal * logMi))
  }
  probs
}

tg_surps <- function(theta, analyze_res) {
  iter <- length(analyze_res$parList)
  no_items <- length(analyze_res$parList[[iter]]$WfdList)
  surps <- vector("numeric", 0)
  for (item in 1:no_items) {
    WListi <- analyze_res$parList[[iter]]$WfdList[[item]]
    Wfdi <- WListi$Wfd
    Mi <- WListi$M
    logMi <- log(Mi)
    surps <- c(surps, eval.surp(theta, Wfdi))
  }
  surps
}

mirt_probs <- function(theta, mirt_mod) {
  c(probtrace(mirt_mod, theta))
}

mirt_surps <- function(theta, mirt_mod) {
  no_items <- length(coef(mirt_mod))-1
  surps <- vector("numeric", 0)
  for (item in 1:no_items) {
    extr <- extract.item(mirt_mod, item)
    probs <- probtrace(extr, theta)
    surps <-  c(surps, -log(probs)/log(ncol(probs)))
  }
  surps
}

