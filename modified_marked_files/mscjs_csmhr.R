mscjs_csmhr = function (x, ddl, dml, model_data = NULL, parameters, accumulate = TRUE, 
                        initial = NULL, method, hessian = FALSE, debug = FALSE, chunk_size = 1e+07, 
                        refit, itnmax = NULL, control = NULL, scale, re = FALSE, 
                        compile = FALSE, extra.args = "", clean = TRUE, ...) 
{
  accumulate = FALSE
  nocc = x$nocc
  if (!is.null(ddl$Phi$time.interval)) 
    time.intervals = matrix(ddl$Phi$time.interval, nrow(x$data), 
                            ncol = nocc - 1, byrow = TRUE)
  else if (is.vector(x$time.intervals)) 
    time.intervals = matrix(x$time.intervals, nrow = nrow(x$data), 
                            ncol = nocc - 1, byrow = TRUE)
  else time.intervals = x$time.intervals
  strata.labels = x$strata.labels
  uS = x$unobserved
  x = x$data
  freq = NULL
  if (!is.null(x$freq)) 
    freq = x$freq
  ch = x$ch
  imat = process.ch(ch, freq, all = FALSE)
  chmat = matrix((unlist(strsplit(ch, ","))), byrow = TRUE, 
                 ncol = nocc, nrow = length(ch))
  for (nlabel in 1:length(strata.labels)) chmat = t(apply(chmat, 
                                                          1, sub, pattern = strata.labels[nlabel], replacement = nlabel))
  if (is.null(initial)) 
    par = list(Psi = rep(0, ncol(dml$Psi$fe)), p = rep(0, 
                                                       ncol(dml$p$fe)), S = rep(0, ncol(dml$S$fe)))
  else par = set.initial(names(dml), dml, initial)$par
  model_data = list(S.dm = dml$S$fe, p.dm = dml$p$fe, Psi.dm = dml$Psi$fe, 
                    imat = imat, S.fixed = parameters$S$fixed, p.fixed = parameters$p$fixed, 
                    Psi.fixed = parameters$Psi$fixed, time.intervals = time.intervals)
  if (accumulate) {
    cat("Accumulating capture histories based on design. This can take awhile.\n")
    flush.console()
    model_data.save = model_data
  }
  else model_data.save = model_data
  scale = 1
  scale = marked:::set.scale(names(dml), model_data, scale)
  model_data = marked:::scale.dm(model_data, scale)
  if (!re) 
    tpl = "multistate_hazrates"
  else stop("random effect portion not completed for this model")
  marked:::setup_admb(tpl, compile, clean, re = FALSE)
  con = file(paste(tpl, ".dat", sep = ""), open = "wt")
  n = length(model_data$imat$freq)
  write(n, con, append = FALSE)                                        # n
  nocc = model_data$imat$nocc
  write(nocc, con, append = TRUE)                                      # m
  nS = length(strata.labels)
  write(nS, con, append = TRUE)                                        # nS
  write(t(chmat), con, ncolumns = nocc, append = TRUE)                 # CH
  write(model_data$imat$first, con, ncolumns = n, append = TRUE)       # first
  if (!re) {
    write(model_data$imat$freq, con, ncolumns = n, append = TRUE)    # freq
  }
  else {
    if (any(model_data$imat$freq != 1)) 
      stop("\n cannot use random effects with frequency >1")
  }
  write(t(model_data$time.intervals), con, ncolumns = nocc - 
          1, append = TRUE)                                                # tint
  phidm = as.matrix(model_data$S.dm)
  phifix = rep(-1, nrow(phidm))
  if (!is.null(ddl$S$fix)) 
    phifix[!is.na(ddl$S$fix)] = ddl$S$fix[!is.na(ddl$S$fix)]
  slist = marked:::simplify_indices(cbind(phidm, phifix))
  write(ncol(phidm), con, append = TRUE)                                      # kphi
  write(length(slist$set), con, append = TRUE)                                # nrowphi
  #    write(t(phidm[slist$set, , drop = FALSE]), con, ncolumns = ncol(phidm), 
  #        append = TRUE)
  write(phifix[slist$set], con, append = TRUE)                                # phifix
  write(slist$indices[ddl$S.indices], con, ncolumns = n, append = TRUE)                     # phiindex(1,all_nrows);      // phi indices
  pdm = as.matrix(model_data$p.dm)
  pfix = rep(-1, nrow(pdm))
  if (!is.null(ddl$p$fix)) 
    pfix[!is.na(ddl$p$fix)] = ddl$p$fix[!is.na(ddl$p$fix)]
  slist = marked:::simplify_indices(cbind(pdm, pfix))
  write(ncol(pdm), con, append = TRUE)                                        # kp
  write(length(slist$set), con, append = TRUE)                                # nrowp
  write(t(pdm[slist$set, , drop = FALSE]), con, ncolumns = ncol(pdm), 
        append = TRUE)                                                          # pdm                                                  
  write(pfix[slist$set], con, append = TRUE)                                  # pfix
  write(slist$indices[ddl$p.indices], con, ncolumns = n, append = TRUE)       # pindex
  psidm = as.matrix(model_data$Psi.dm)
  psifix = rep(-1, nrow(psidm))
  if (!is.null(ddl$Psi$fix)) 
    psifix[!is.na(ddl$Psi$fix)] = ddl$Psi$fix[!is.na(ddl$Psi$fix)]
  slist = marked:::simplify_indices(cbind(psidm, psifix))
  write(ncol(psidm), con, append = TRUE)                                      # kpsi
  write(length(slist$set), con, append = TRUE)                                # nrowpsi
  write(t(psidm[slist$set, , drop = FALSE]), con, ncolumns = ncol(psidm), 
        append = TRUE)
  write(psifix[slist$set], con, append = TRUE)
  write(slist$indices[ddl$Psi.indices], con, ncolumns = n, append = TRUE)
  close(con)
  
  con = file(paste(tpl, ".pin", sep = ""), open = "wt")
  #    write(par$S, con, ncolumns = length(par$S), append = FALSE)
  write(par$p, con, ncolumns = length(par$p), append = TRUE)
  write(par$Psi, con, ncolumns = length(par$Psi), append = FALSE)
  close(con)
  if (hessian) 
    xx = run_admb(tpl, extra.args = extra.args)
  else xx = run_admb(tpl, extra.args = paste(extra.args, "-nohess"))
  convergence = attr(xx, "status")
  if (is.null(convergence)) 
    convergence = 0
  res = read_admb(tpl)
  beta = list(marked:::unscale.par(c(res$coeflist$phibeta, res$coeflist$pbeta, 
                                     res$coeflist$psibeta), scale))
  parnames = names(unlist(beta))
  fixed.npar = length(unlist(beta))
  if (!is.null(res$hes)) {
    beta.vcv = marked:::solvecov(res$hes)$inv
    rownames(res$hes) = parnames
    colnames(res$hes) = rownames(res$hes)
    if (all(diag(beta.vcv > 0))) 
      res$cor = beta.vcv/outer(sqrt(diag(beta.vcv)), sqrt(diag(beta.vcv)))
  }
  else beta.vcv = res$vcov
  rownames(beta.vcv) = parnames
  colnames(beta.vcv) = rownames(beta.vcv)
  rownames(res$cor) = rownames(beta.vcv)
  colnames(res$cor) = rownames(beta.vcv)
  res$vcov = NULL
  optim.details = c(fn = res$fn, maxgrad = res$maxgrad, eratio = res$eratio)
  options = list(extra.args = extra.args)
  res$cor = NULL
  res$maxgrad = NULL
  results = c(beta = beta, neg2lnl = -2 * res$loglik, AIC = -2 * 
                res$loglik + 2 * res$npar, convergence = convergence)
  results$optim.details = optim.details
  results$options = options
  results$coeflist = res$coeflist
  results$npar = list(npar = res$npar, npar_sdrpt = res$npar_sdrpt, 
                      npar_total = res$npar_total)
  results$beta.vcv = beta.vcv
  res = results
  res$model_data = model_data.save
  class(res) = c("crm", "admb", "mle", "mscjs")
  return(res)
}