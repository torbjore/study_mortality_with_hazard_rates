crm_csmhr = function (data, ddl = NULL, begin.time = 1, model = "MSCJS", title = "", 
                      model.parameters = list(), design.parameters = list(), initial = NULL, 
                      groups = NULL, time.intervals = NULL, debug = FALSE, method = "BFGS", 
                      hessian = FALSE, accumulate = TRUE, chunk_size = 1e+07, control = list(), 
                      refit = 1, itnmax = 5000, scale = NULL, run = TRUE, burnin = 100, 
                      iter = 1000, use.admb = FALSE, use.tmb = FALSE, crossed = NULL, 
                      reml = FALSE, compile = FALSE, extra.args = NULL, strata.labels = NULL, 
                      clean = NULL, save.matrices = TRUE, simplify = FALSE, getreals = FALSE, 
                      check = FALSE, csmhr=TRUE, ...) 
{
  model = toupper(model)
  ptm = proc.time()
  if (is.null(crossed)) 
    crossed = FALSE
  if (crossed) 
    accumulate = FALSE
  if (is.null(data$data)) {
    if (!is.null(ddl)) {
      warning("Warning: specification of ddl ignored, as data have not been processed")
      ddl = NULL
    }
    message("Model: ", model, "\n")
    message("Processing data...\n")
    flush.console()
    data.proc = process.data(data, begin.time = begin.time, 
                             model = model, mixtures = 1, groups = groups, age.var = NULL, 
                             initial.ages = NULL, time.intervals = time.intervals, 
                             nocc = NULL, accumulate = accumulate, strata.labels = strata.labels)
  }
  else {
    data.proc = data
    model = data$model
  }
  if(csmhr)
  {
    states = levels(ddl$Psi$stratum)
    fst = states[1]
    lst = states[length(states)]
    message("Alive state: ", fst, "; Newly dead states: ", states[2:(length(states)-1)])
    ddl$Psi$fix[ddl$Psi$stratum == fst & ddl$Psi$tostratum == lst] = 0
    ddl$Psi$fix[ddl$Psi$stratum != fst] = 1 # NOTE: It does not matter what values you fix to here as the value is reset in the tpl-file by mscjs_csmhr
    ddl$p$fix[ddl$p$stratum==lst] = 0
    ddl$S$fix[ddl$S$stratum!=lst] = 1
    ddl$S$fix[ddl$S$stratum==lst] = 0
    if(!is.null(model.parameters[["S"]])) warning("Model specification for S is ignored when csmhr == TRUE. Cause specific mortality hazard rates are instead modelled as transition parameters from the 1st state to the other states (Psi = log hazard rates)")
    model.parameters[["S"]] = list(formula = ~ 0)
  }
  number.of.groups = 1
  if (!is.null(data.proc$group.covariates)) 
    number.of.groups = nrow(data.proc$group.covariates)
  par.list = setup.parameters(data.proc$model, check = TRUE)
  if (!marked:::valid.parameters(model, model.parameters)) 
    stop()
  parameters = setup.parameters(data.proc$model, model.parameters, 
                                data$nocc, number.of.groups = number.of.groups)
  parameters = parameters[par.list]
  re = FALSE
  for (i in 1:length(parameters)) {
    if (is.null(parameters[[i]]$formula)) 
      parameters[[i]]$formula = ~1
    mlist = marked:::proc.form(parameters[[i]]$formula)
    if (!is.null(mlist$re.model)) {
      re_names = sub("^\\s+", "", sapply(strsplit(names(mlist$re.model), 
                                                  "\\|"), function(x) x[2]))
      if (length(re_names) > 1 | !"id" %in% re_names) 
        crossed = TRUE
      if ((length(re_names) > 1 || re_names[1] != "time" || 
           use.admb) & any(data.proc$freq > 1)) 
        stop("\n data cannot be accumulated (freq>1) except with temporal random effects only; set accumulate=FALSE\n")
      re = TRUE
    }
    if (parameters[[i]]$nointercept) 
      parameters[[i]]$remove.intercept = TRUE
  }
  if (re & !use.tmb) {
    use.admb = TRUE
    if (is.null(clean)) 
      clean = TRUE
  }
  if (use.admb) {
    if (!re) 
      crossed = FALSE
    if (is.null(clean)) 
      clean = TRUE
  }
  if (use.tmb & is.null(clean)) 
    clean = FALSE
  if (is.null(ddl)) {
    message("Creating design data...\n")
    flush.console()
    ddl = make.design.data(data.proc, design.parameters)
  }
  else {
    for (i in 1:length(parameters)) {
      if (!is.null(ddl[[i]]$order)) 
        if (any(ddl[[i]]$order != 1:nrow(ddl[[i]]))) 
          stop(paste("Design data for parameter", names(parameters)[i], 
                     "is out of order."))
    }
    if (!is.null(design.parameters)) 
      for (parname in names(design.parameters)) ddl$design.parameters[[parname]] = c(ddl$design.parameters[[parname]], 
                                                                                     design.parameters[[parname]])
      design.parameters = ddl$design.parameters
  }
  ddl = marked:::set.fixed(ddl, parameters)
  if (model == "CSMHR" | model == "MSCJS" | (substr(model, 1, 4) == "MVMS" & use.admb)) 
    ddl = marked:::simplify_ddl(ddl, parameters)
  if (substr(model, 1, 4) == "MVMS") {
    if (is.null(ddl$pi$fix)) 
      message("\n No values provided for fix for pi. At least need to set a reference cell")
    else {
      bad_pi = sapply(split(ddl$pi$fix, ddl$pi$id), function(x) {
        ifelse(any(is.na(x)) & !(any(x[!is.na(x)] == 
                                       1)), TRUE, FALSE)
      })
      if (any(bad_pi)) 
        message("\n Check values of fix for pi. Reference cell (fix=1) should be set if any are estimated (fix=NA)")
    }
    if (is.null(ddl$delta$fix)) {
      message("\n No values provided for fix for delta. At least need to set a reference cell")
    }
    else {
      bad_delta = sapply(split(ddl$delta$fix, list(ddl$delta$id, 
                                                   ddl$delta$occ, ddl$delta$stratum)), function(x) {
                                                     ifelse(any(is.na(x)) & !(any(x[!is.na(x)] == 
                                                                                    1)), TRUE, FALSE)
                                                   })
      if (any(bad_delta)) 
        message("\n Check values of fix for delta. Reference cell (fix=1) should be set if any are estimated (fix=NA)")
    }
    if (is.null(ddl$Psi$fix)) {
      message("\n No values provided for fix for delta. At least need to set a reference cell")
    }
    else {
      bad_Psi = sapply(split(ddl$Psi$fix, list(ddl$Psi$id, 
                                               ddl$Psi$occ, ddl$Psi$stratum)), function(x) {
                                                 ifelse(any(is.na(x)) & !(any(x[!is.na(x)] == 
                                                                                1)), TRUE, FALSE)
                                               })
      if (any(bad_Psi)) 
        message("\n Check values of fix for Psi. Reference cell (fix=1) should be set if any are estimated (fix=NA)")
    }
  }
  if (simplify) {
    simplify = FALSE
    message("simplify argument has been disabled")
  }
  for (i in 1:length(parameters)) {
    if (!is.null(ddl[[i]]$fix)) {
      if (all(!is.na(ddl[[i]]$fix))) {
        message(paste("All values for", names(parameters)[i], 
                      "have been fixed. Setting formula to ~0"))
        parameters[[i]]$formula = ~0
      }
      else {
        if (parameters[[i]]$formula == ~0) 
          stop(paste("Cannot use formula ~0 for", names(parameters)[i], 
                     "when some of the parameters must be estimated"))
      }
    }
    else if (parameters[[i]]$formula == ~0) 
      stop(paste("Cannot use formula ~0 for", names(parameters)[i], 
                 "when some of the parameters must be estimated"))
  }
  dml = marked:::create.dml(ddl, model.parameters = parameters, design.parameters = design.parameters, 
                            chunk_size = chunk_size, simplify = simplify, use.admb = use.admb)
  if (substr(model, 1, 3) == "HMM" | (nchar(model) >= 4 & substr(model, 
                                                                 1, 4) == "MVMS")) 
    initial.list = set.initial(names(dml), dml, initial)
  else initial.list = NULL
  if (!run) 
    return(list(model = model, data = data.proc, model.parameters = parameters, 
                design.parameters = design.parameters, ddl = ddl, 
                dml = dml, results = initial.list))
  if ("SANN" %in% method) {
    if (length(method) > 1) 
      warning("***SANN can only be used by itself; other methods ignored.")
    method = "SANN"
    control$maxit = itnmax
  }
  if ("nlminb" %in% method) {
    control$eval.max = itnmax
    control$iter.max = itnmax
  }
  message("Fitting model\n")
  if (model == "CJS") 
    if (use.tmb) {
      runmodel = cjs_tmb(data.proc, ddl, dml, parameters = parameters, 
                         initial = initial, method = method, hessian = hessian, 
                         debug = debug, accumulate = accumulate, chunk_size = chunk_size, 
                         refit = refit, control = control, itnmax = itnmax, 
                         scale = scale, crossed = crossed, compile = compile, 
                         extra.args = extra.args, reml = reml, clean = clean, 
                         getreals = getreals, ...)
    }
  else runmodel = cjs(data.proc, ddl, dml, parameters = parameters, 
                      initial = initial, method = method, hessian = hessian, 
                      debug = debug, accumulate = accumulate, chunk_size = chunk_size, 
                      refit = refit, control = control, itnmax = itnmax, 
                      scale = scale, use.admb = use.admb, crossed = crossed, 
                      compile = compile, extra.args = extra.args, reml = reml, 
                      clean = clean, ...)
  if (model == "JS") 
    runmodel = js(data.proc, ddl, dml, parameters = parameters, 
                  initial = initial, method = method, hessian = hessian, 
                  debug = debug, accumulate = FALSE, chunk_size = chunk_size, 
                  refit = refit, control = control, itnmax = itnmax, 
                  scale = scale, ...)
  if (model == "MSCJS" & !csmhr) 
    runmodel = mscjs(data.proc, ddl, dml, parameters = parameters, 
                     initial = initial, method = method, hessian = hessian, 
                     debug = debug, accumulate = accumulate, chunk_size = chunk_size, 
                     refit = refit, control = control, itnmax = itnmax, 
                     scale = scale, re = re, compile = compile, extra.args = extra.args, 
                     clean = clean, ...)
  if (model == "MSCJS" & csmhr) 
    runmodel = mscjs_csmhr(data.proc, ddl, dml, parameters = parameters, 
                           initial = initial, method = method, hessian = hessian, 
                           debug = debug, accumulate = accumulate, chunk_size = chunk_size, 
                           refit = refit, control = control, itnmax = itnmax, 
                           scale = scale, re = re, compile = compile, extra.args = extra.args, 
                           clean = clean, ...)
  if (model == "PROBITCJS") {
    if (is.null(initial)) {
      imat = process.ch(data.proc$data$ch, data.proc$data$freq, 
                        all = FALSE)
      runmodel = probitCJS(ddl, dml, parameters = parameters, 
                           design.parameters = design.parameters, imat = imat, 
                           iter = iter, burnin = burnin)
    }
    else runmodel = probitCJS(ddl, dml, parameters = parameters, 
                              design.parameters = design.parameters, initial = initial, 
                              iter = iter, burnin = burnin)
  }
  if (substr(model, 1, 3) == "HMM" | (nchar(model) >= 4 & substr(model, 
                                                                 1, 4) == "MVMS")) {
    if (substr(model, 1, 4) == "MVMS") 
      sup = data.proc$fct_sup(list(obslevels = data.proc$ObsLevels))
    else sup = NULL
    if (is.null(data.proc$strata.list) | substr(model, 1, 
                                                4) == "MVMS") {
      mx = data.proc$m
    }
    else {
      mx = list(ns = length(data.proc$strata.list$states), 
                na = length(data.proc$strata.list[[names(data.proc$strata.list)[names(data.proc$strata.list) != 
                                                                                  "states"]]]))
    }
    if (use.admb & model == "MVMSCJS") {
      xx = HMMLikelihood(par = unlist(initial.list$par), 
                         xx = data.proc$ehmat, mx = mx, type = initial.list$ptype, 
                         T = data.proc$nocc, xstart = data.proc$start, 
                         freq = data.proc$freq, fct_dmat = data.proc$fct_dmat, 
                         fct_gamma = data.proc$fct_gamma, fct_delta = data.proc$fct_delta, 
                         ddl = ddl, dml = dml, parameters = parameters, 
                         sup = sup, check = TRUE)
      runmodel = mvmscjs(data.proc, ddl, dml, parameters = parameters, 
                         initial = initial, method = method, hessian = hessian, 
                         debug = debug, accumulate = accumulate, chunk_size = chunk_size, 
                         refit = refit, control = control, itnmax = itnmax, 
                         scale = scale, re = re, compile = compile, extra.args = extra.args, 
                         clean = clean, sup = sup, ...)
      par = coef(runmodel)[, 1]
      runmodel$options = c(runmodel$options, list(accumulate = accumulate, 
                                                  initial = initial.list$par, method = method, 
                                                  chunk_size = chunk_size, itnmax = itnmax, control = control))
    }
    else {
      xx = HMMLikelihood(par = unlist(initial.list$par), 
                         xx = data.proc$ehmat, mx = mx, type = initial.list$ptype, 
                         T = data.proc$nocc, xstart = data.proc$start, 
                         freq = data.proc$freq, fct_dmat = data.proc$fct_dmat, 
                         fct_gamma = data.proc$fct_gamma, fct_delta = data.proc$fct_delta, 
                         ddl = ddl, dml = dml, parameters = parameters, 
                         sup = sup, check = TRUE)
      runmodel = optimx(unlist(initial.list$par), HMMLikelihood, 
                        method = method, debug = debug, hessian = hessian, 
                        itnmax = itnmax, xx = data.proc$ehmat, mx = mx, 
                        type = initial.list$ptype, T = data.proc$nocc, 
                        xstart = data.proc$start, freq = data.proc$freq, 
                        control = control, fct_dmat = data.proc$fct_dmat, 
                        fct_gamma = data.proc$fct_gamma, fct_delta = data.proc$fct_delta, 
                        ddl = ddl, dml = dml, parameters = parameters, 
                        sup = sup, check = FALSE)
      par = coef(runmodel, order = "value")[1, ]
      runmodel = list(optim.details = as.list(summary(runmodel, 
                                                      order = "value", par.select = FALSE)[1, ]))
      if (hessian) 
        runmodel$hessian = attr(runmodel$optim.details, 
                                "details")$nhatend
      runmodel$convergence = runmodel$optim.details$convcode
      runmodel$options = list(accumulate = accumulate, 
                              initial = initial.list$par, method = method, 
                              chunk_size = chunk_size, itnmax = itnmax, control = control)
    }
    if (save.matrices) {
      runmodel$mat = HMMLikelihood(par = par, type = initial.list$ptype, 
                                   xx = data.proc$ehmat, mx = mx, T = data.proc$nocc, 
                                   xstart = data.proc$start, freq = data.proc$freq, 
                                   fct_dmat = data.proc$fct_dmat, fct_gamma = data.proc$fct_gamma, 
                                   fct_delta = data.proc$fct_delta, ddl = ddl, dml = dml, 
                                   parameters = parameters, return.mat = TRUE, sup = sup)
      if (model == "HMMCJS") {
        dimnames(runmodel$mat$gamma)[3:4] = list(c("Alive", 
                                                   "Dead"), c("Alive", "Dead"))
        dimnames(runmodel$mat$dmat)[3:4] = list(c("Missed", 
                                                  "Seen"), c("Alive", "Dead"))
      }
      else {
        dimnames(runmodel$mat$gamma)[3:4] = list(c(data.proc$strata.labels, 
                                                   "Dead"), c(data.proc$strata.labels, "Dead"))
        dimnames(runmodel$mat$dmat)[3:4] = list(data.proc$ObsLevels, 
                                                c(data.proc$strata.labels, "Dead"))
      }
      names(dimnames(runmodel$mat$gamma)) = c("Id", "Occasion", 
                                              "From_state", "To_state")
      names(dimnames(runmodel$mat$dmat)) = c("Id", "Occasion", 
                                             "Observation", "State")
    }
    parlist = split(par, initial.list$ptype)
    par = vector("list", length = length(names(initial.list$par)))
    names(par) = names(initial.list$par)
    for (p in names(parlist)) {
      par[[p]] = parlist[[p]]
      names(par[[p]]) = colnames(dml[[p]]$fe)
    }
    runmodel$beta = par
    runmodel$par = NULL
    if (is.null(runmodel$neg2lnl)) 
      runmodel$neg2lnl = 2 * runmodel$optim.details$value
    runmodel$AIC = runmodel$neg2lnl + 2 * sum(sapply(runmodel$beta, 
                                                     length))
    if (!is.null(runmodel$hessian)) {
      runmodel$beta.vcv = solvecov(runmodel$hessian)$inv
      colnames(runmodel$beta.vcv) = names(unlist(runmodel$beta))
      rownames(runmodel$beta.vcv) = colnames(runmodel$beta.vcv)
    }
    class(runmodel) = c("crm", "mle", model)
  }
  if (!is.null(runmodel$convergence) && runmodel$convergence != 
      0 & !use.admb) {
    warning("******Model did not converge******")
    msg = attr(runmodel$optim.details, "details")$message
    if (is.null(msg)) 
      msg = "Exceeded maximum number of iterations"
    warning(msg)
  }
  object = list(model = model, data = data.proc, model.parameters = parameters, 
                design.parameters = design.parameters, results = runmodel)
  class(object) = class(runmodel)
  if (!re & model != "MSCJS" & (nchar(model) < 4 | (nchar(model) >= 
                                                    4 & substr(model, 1, 4) != "MVMS"))) 
    object$results$reals = predict(object, ddl = ddl, unique = TRUE, 
                                   se = hessian)
  cat(paste("\nElapsed time in minutes: ", round((proc.time()[3] - 
                                                    ptm[3])/60, digits = 4), "\n"))
  return(object)
}
