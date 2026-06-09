###############
## Libraries
###############
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(doMC)))
suppressWarnings(suppressMessages(library(multcomp)))
suppressWarnings(suppressMessages(library(argparser)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(lme4)))

# IMPROVED: Add better convergence control
glmerControl_settings <- glmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 100000),
  calc.derivs = FALSE,
  tolPwrss = 1e-3
)

# Alternative optimizers to try if default fails
alt_optimizers <- list(
  bobyqa = list(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
  nloptwrap = list(optimizer = "nloptwrap", 
                   optCtrl = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 100000)),
  Nelder_Mead = list(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 100000))
)

## FUNCTIONS
#################
parse_cl_args=function(args){
  # Create a parser
  p <- arg_parser("Calculate GLM coefficients and p-values")
  
  p <- add_argument(p, "HAFs", help=".Rdata file containing 3 objects: afmat(matrix), samples(data.frame), sites(data.frame)")
  p <- add_argument(p, "--readDepth", help=".RDS file containing a matrix with same dimensions as afmat in HAFs, giving the raw read depth per site/sample", default=NA)
  p <- add_argument(p, "--effectiveCov", help="either a single number to be used as the effective coverage for every site/sample, or \na .csv file with column names ['sampID','chrom','ec'] containing an estimate of effective coverage per chrom/sample", default=NA)
  p <- add_argument(p, "--chrom", help="run GLM only for sites on this chromosome", default=NA)
  p <- add_argument(p, "--dropRep", help="ID of replicate to drop (when running leave-one-out)", default=NA)
  p <- add_argument(p, "--poolSize", help="number individuals sampled per pool", default=100,type="integer")
  p <- add_argument(p, "--mainEffect", help="calculate p-values for all pairwise comparisons of groups in this sample metadata column")
  p <- add_argument(p, "--effect1", help="calculate p-values for a variable1 without pairwise comparisons")
  p <- add_argument(p, "--effect1type", help="type of data (continuous, categorical)")  
  p <- add_argument(p, "--effect2", help="calculate p-values for a variable2 without pairwise comparisons")
  p <- add_argument(p, "--effect2type", help="type of data (continuous, categorical)")
  p <- add_argument(p, "--interaction", help="calculate p-values for an interaction of effect1 and effect2. Options: TRUE or FALSE")
  p <- add_argument(p, "--repName", help="name of the column in the sample metadata table that identifies replicate IDs")
  p <- add_argument(p, "--repNamecontrast", help="sum, alt")
  p <- add_argument(p, "--repNameignore", help="add rep to the linear model or not? TRUE or FALSE")
  p <- add_argument(p, "--repNameinter", help="add an interaction with rep? If so, with which effect? options: TRUE or FALSE")
  p <- add_argument(p, "--repNameeff", help="which effect to use as an interactor with rep? options: effect1 or effect2")
  p <- add_argument(p, "--testNsites", help="run GLM on a random subset of N sites", default=NA,type="integer")
  p <- add_argument(p, "--nCores", help="run GLM in parallel using mclapply with this many cores", default=2,type="integer")
  p <- add_argument(p, "--saveAs", help="format for saving results dataframe: 'RDS', 'Rdata', or 'csv'", default="RDS")
  p <- add_argument(p, "--outDir", help="write all results to this directory; will be created if it doesnt already exist", default=".")
  p <- add_argument(p, "--trtContrast", help="specify reference level for treatment contrast. Example: M", default=NA)
  
  p <- add_argument(p, "--selectTrts", help="Comma-separated list of treatments to keep (e.g. 'B,T'). If NA, uses all.", default=NA)
  
  p <- add_argument(p, "--makePlots", help="create plots for each GLM: TRUE or FALSE", default="FALSE")  
  p <- add_argument(p, "--saveModels", help="save GLM model objects: TRUE or FALSE", default="FALSE") 
  p <- add_argument(p, "--trtPopGroups", help="create treatment-population groups (trt_pop): TRUE or FALSE", default="FALSE")
  p <- add_argument(p, "--mixedEffects", help="use mixed-effects model with replicate random effects: TRUE or FALSE", default="FALSE")
  p <- add_argument(p, "--compareModels", help="compare multiple mixed-effects models: TRUE or FALSE", default="FALSE")
  p <- add_argument(p, "--randomEffectsType", help="type of random effects: simple, nested, crossed", default="simple")

  args <- parse_args(p)
  
  ## check if the HAF files exist
  if(! file.exists(args$HAFs) ){cat("HAFs file",args$HAFs,"does not exist\n***EXITING***\n");quit() }
  if(! grepl("\\.Rdata$",args$HAFs) ){cat("HAFs file",args$HAFs,"must be saved as an .Rdata file\n***EXITING***\n");quit() }
  
  ## must supply either readDepth or effectiveCov; effCov can be a file or a single integer
  ## readDepth should be a file containing a matrix with the same dimensions as afmat in HAFs, giving the raw read depth per site/sample
  ## if effectiveCov is a number, it will be used as the effective coverage for every site/sample
  ## if effectiveCov is a file, it should contain an estimate of effective coverage per chrom/sample

  if(is.na(args$readDepth)){
    if(is.na(args$effectiveCov)){
      cat("either --readDepth or --effectiveCov must be supplied.\n***EXITING***\n");quit()
    } else{
      if(is.na(as.numeric(args$effectiveCov)) & !file.exists(args$effectiveCov)){
        cat("--effectiveCov must either be a single integer or a file\n***EXITING***\n");quit()
      } else {
        if(!file.exists(args$effectiveCov)){args$effectiveCov=as.integer(args$effectiveCov)}
      }
    }
  } else {
    if(!file.exists(args$readDepth)){cat("readDepth file",args$readDepth,"does not exist\n***EXITING***\n");quit()}
  }
  
  ## check nCores,testNsites,poolSize
  if (is.na(as.numeric(args$nCores))) {cat("nCores must be an integer\n***EXITING***\n");quit()}
  if (!is.na(args$testNsites) && is.na(as.numeric(args$testNsites))) {cat("testNsites must be an integer\n***EXITING***\n");quit()}
  if (is.na(args$poolSize)) {cat("poolSize must be an integer\n***EXITING***\n");quit()}
 
  if (args$mainEffect == "NA") {args$mainEffect=NULL}
  if (args$effect1 == "NA") {args$effect1=NA}
  if (args$effect1type == "NA") {args$effect1type=NA}
  if (args$effect2 == "NA") {args$effect2=NA}
  if (args$effect2type == "NA") {args$effect2type=NA}
  if (args$repNamecontrast == "NA") {args$repNamecontrast=NA}

  return(args)
}

check_convergence <- function(model) {
  warnings <- model@optinfo$conv$lme4$messages
  singular <- isSingular(model)
  return(list(converged = is.null(warnings), singular = singular))
}

check_overdispersion <- function(model) {
  pearson_resid <- residuals(model, type = "pearson")
  phi <- sum(pearson_resid^2) / df.residual(model)
  return(phi)
}

set_up_sampData=function(sampIX, samps, model.vars, repName, repNamecontrast, cmpAll, effect1, effect1type, effect2, effect2type, trtContrast, trtPopGroups=FALSE){
  sampData=samps[sampIX,]%>%dplyr::select(all_of(model.vars))
  
  # Ensure factors are clean (drop unused levels from removed treatments)
  sampData[] <- lapply(sampData, function(x) if(is.factor(x)) droplevels(x) else x)

  varCols=which(apply(sampData,2,function(x){length(x)>1}))
  sampData=sampData%>%dplyr::select(all_of(varCols))
  
  print(varCols)
  
  if(trtPopGroups && "trt" %in% names(sampData) && repName %in% names(sampData)){
    sampData$trt_pop <- interaction(sampData$trt, sampData[[repName]], sep="_", drop=TRUE)
    sampData$trt_pop <- factor(sampData$trt_pop)
    cat("Created trt_pop groups:", levels(sampData$trt_pop), "\n")
  }
  
  if(!is.na(trtContrast) && "trt" %in% names(sampData)){
    sampData$trt = factor(sampData$trt)
    all_levels = levels(sampData$trt)
    new_levels = c(trtContrast, all_levels[all_levels != trtContrast])
    sampData$trt = factor(sampData$trt, levels = new_levels)
    cat("Treatment contrasts set with", trtContrast, "as reference level\n")
  }
  
  if(!is.null(cmpAll)){
    sampData[[cmpAll]]=factor(sampData[[cmpAll]])
  }
    
  if(!is.na(effect1)){
    if(effect1type == "continuous"){
        sampData[, effect1] = as.numeric(sampData[, effect1])
    } else if(effect1type == "categorical"){
        sampData[, effect1] = as.factor(sampData[, effect1])
    } else { print("Need to specify data type")}
  }
    
  if(!is.na(effect2)){
    if(effect2type == "continuous"){
        sampData[, effect2] = as.numeric(sampData[, effect2])
    } else if(effect2type == "categorical"){
        sampData[, effect2] = as.factor(sampData[, effect2])
    } else { print("Need to specify data type")}
  }

  if(!is.na(repName)){
    if(repNamecontrast == "sum"){
        sampData[, repName] = as.factor(sampData[, repName])
        contrasts(sampData[, repName]) <- matrix(c(
          1,  0,  0,
         -1,  1,  0,
          0, -1,  1,
          0,  0, -1
        ), nrow = 4, byrow = TRUE)
    } else if(repNamecontrast == "alt"){
        # TBD
    } else { 
        sampData[, repName] = as.factor(sampData[, repName])
    }
  }
 
  return(sampData)
} 

set_up_depthData=function(args,sites,samps){
  if(!is.na(args$readDepth)){
    cat("using raw read depth for binomial counts\n")
    rd=data.frame(readRDS(args$readDepth))
    ss=match(samps$SourceFile,colnames(rd))
    if(sum(is.na(ss)>0)){cat("--readDepth file is missing samples:",unique(samps$SourceFile[is.na(ss)]),"\n***EXITING***\n");quit()}
    rd=rd[,ss]
    if(nrow(data.frame(rd))!=nrow(data.frame(sites))){cat("--readDepth file must have the same number of rows as sites/afmat in --HAFs file\n***EXITING***\n");quit()}
  } else{
    if(!is.numeric(args$effectiveCov)){
      cat("using per-chromosome effective coverage for binomial counts\n")
      df.ec=fread(args$effectiveCov);
      if(length(intersect(c("sampID","chrom","ec"),colnames(df.ec)))<3){
        cat("--effectiveCov file must have the following column names: sampID,chrom,ec\n***EXITING***\n");quit()
      }
      df.ec$chrom <- factor(df.ec$chrom)
      mat.ec = t(xtabs(ec ~ sampID + chrom, df.ec))
      cc=match(sites$chrom,rownames(mat.ec))
      ss=match(samps$SourceFile,colnames(mat.ec))
      if(sum(is.na(cc))>0){cat("--effectiveCov file is missing chromosomes:",unique(sites$chrom[is.na(cc)]),"\n***EXITING***\n");quit()}
      if(sum(is.na(ss))>0){cat("--effectiveCov file is missing samples:",unique(samps$SourceFile[is.na(ss)]),"\n***EXITING***\n");quit()}
      rd=round(mat.ec[cc,ss])
    } else{
      cat("using constant value for binomial counts\n")
      rd=matrix(args$effectiveCov,nrow(sites),nrow(samps))
    }
  }
  return(rd)
}

extract_coef_pval=function(model,cmpAll=NULL,dontReport=NULL,includeIntercept=TRUE){
  if(!is.null(cmpAll)){
    model.multcomp=summary(eval(parse(text=paste0("glht(model, mcp(",cmpAll,"='Tukey'))")))) 
    cp=cbind(coefficients(model.multcomp),model.multcomp$test$pvalues)
    row.names(cp)=gsub(" - ","_",row.names(cp))
  } else{
    if(includeIntercept) {
      cp=summary(model)$coefficients[,c(1,4),drop=FALSE];
    } else {
    cp=summary(model)$coefficients[-1,c(1,4),drop=FALSE];
    cp=summary(model)$coefficients[-1,c(1,4),drop=FALSE];
    }
  }
  if(!is.null(dontReport)){
    cp=cp[grep(paste0("(",paste0(dontReport,collapse="|"),")"),row.names(cp),invert = TRUE),,drop=FALSE];
  }
  
  return(cp)
}

fit_glmer_robust <- function(formula, data, family = binomial) {
  model <- NULL
  last_error <- NULL
  
  for(opt_name in names(alt_optimizers)) {
    try({
      cat("Trying optimizer:", opt_name, "\n")
      ctrl <- do.call(glmerControl, alt_optimizers[[opt_name]])
      model <- glmer(formula, family = family, data = data, control = ctrl)
      
      conv_check <- check_convergence(model)
      if(conv_check$converged && !conv_check$singular) {
        cat("Success with optimizer:", opt_name, "\n")
        return(model)
      } else if(!is.null(model)) {
        cat("Optimizer", opt_name, "converged but singular\n")
      }
    }, silent = FALSE)
  }
  
  if(!is.null(model)) {
    cat("Returning best available model (may be singular)\n")
    return(model)
  } else {
    cat("All optimizers failed\n")
    return(NULL)
  }
}

fit_multiple_models <- function(af.site, rd.site, sampData, formulaBase, repName, randomEffectsType = "simple") {
  Neff = rd.site
  cts = cbind(round(as.numeric(Neff)*as.numeric(af.site)), 
              round(as.numeric(Neff)*(1-as.numeric(af.site))))
  
  if(!"pop" %in% names(sampData) && repName %in% names(sampData)) {
    sampData$pop <- sampData[[repName]]
  }
  
  if(!"trt_pop" %in% names(sampData) && "trt" %in% names(sampData) && repName %in% names(sampData)) {
    sampData$trt_pop <- interaction(sampData$trt, sampData[[repName]], sep="_", drop=TRUE)
  }
  
  default_result <- data.frame(
    model = NA,
    AIC = NA,
    BIC = NA,
    logLik = NA,
    converged = FALSE,
    singular = TRUE,
    overdispersion = NA,
    formula_used = NA,
    stringsAsFactors = FALSE
  )
  
  results <- list()
  models <- list()
  
  if(randomEffectsType == "simple") {
     formulas <- list(
       "intercepts_only" = "cts ~ gen * trt + (1 | pop)"
  )
  } else if(randomEffectsType == "nested") {
     formulas <- list(
       "intercepts_only" = "cts ~ gen * trt + (1 | trt_pop)"
  )
  } else { # crossed
     formulas <- list(
       "intercepts_only" = "cts ~ gen * trt + (1 | pop) + (1 | trt)"
  )
}

  # Fit each model
  for(model_name in names(formulas)) {
    result <- default_result
    result$model <- model_name
    result$formula_used <- formulas[[model_name]]
    
    model <- fit_glmer_robust(as.formula(formulas[[model_name]]), sampData, binomial)
    
    if(!is.null(model)) {
      models[[model_name]] <- model
      
      conv_check <- check_convergence(model)
      overdispersion <- check_overdispersion(model)
      
      result <- data.frame(
        model = model_name,
        AIC = AIC(model),
        BIC = BIC(model),
        logLik = as.numeric(logLik(model)),
        converged = conv_check$converged,
        singular = conv_check$singular,
        overdispersion = overdispersion,
        formula_used = formulas[[model_name]],
        stringsAsFactors = FALSE
      )
      
      # Add fixed effects
      fixed_effects <- fixef(model)
      fixed_se <- sqrt(diag(vcov(model)))
      fixed_z <- fixed_effects / fixed_se
      fixed_p <- 2 * pnorm(abs(fixed_z), lower.tail = FALSE)
      
      for(i in 1:length(fixed_effects)) {
        coef_name <- paste0("coef_", model_name, ".", make.names(names(fixed_effects)[i]))
        p_name <- paste0("p_", model_name, ".", make.names(names(fixed_effects)[i]))
        result[[coef_name]] <- fixed_effects[i]
        result[[p_name]] <- fixed_p[i]
      }
    }
    
    results[[model_name]] <- result
  }
  
  # Model comparisons (only if models exist)
  comparison_results <- data.frame(
    best_model_by_AIC = NA,
    best_model_by_BIC = NA,
    int_vs_both_chisq = NA,
    int_vs_both_df = NA,
    int_vs_both_p = NA,
    slope_vs_both_chisq = NA,
    slope_vs_both_df = NA,
    slope_vs_both_p = NA,
    stringsAsFactors = FALSE
  )
  
  # Find best models
  valid_models <- results[sapply(results, function(x) !is.na(x$AIC))]
  if(length(valid_models) > 0) {
    aics <- sapply(valid_models, function(x) x$AIC)
    bics <- sapply(valid_models, function(x) x$BIC)
    comparison_results$best_model_by_AIC <- names(valid_models)[which.min(aics)]
    comparison_results$best_model_by_BIC <- names(valid_models)[which.min(bics)]
  }
  
  # Model comparisons
  if("intercepts_only" %in% names(models) && "both" %in% names(models)) {
    try({
      comp1 <- anova(models$intercepts_only, models$both)
      if(nrow(comp1) >= 2) {
        comparison_results$int_vs_both_chisq <- comp1$Chisq[2]
        comparison_results$int_vs_both_df <- comp1$Df[2]
        comparison_results$int_vs_both_p <- comp1$`Pr(>Chisq)`[2]
      }
    }, silent = TRUE)
  }
  
  if("slopes_only" %in% names(models) && "both" %in% names(models)) {
    try({
      comp2 <- anova(models$slopes_only, models$both)
      if(nrow(comp2) >= 2) {
        comparison_results$slope_vs_both_chisq <- comp2$Chisq[2]
        comparison_results$slope_vs_both_df <- comp2$Df[2]
        comparison_results$slope_vs_both_p <- comp2$`Pr(>Chisq)`[2]
      }
    }, silent = TRUE)
  }
  
  # Align all results
  all_cols <- unique(unlist(lapply(results, names)))
  for(i in 1:length(results)) {
    missing_cols <- setdiff(all_cols, names(results[[i]]))
    for(col in missing_cols) {
      results[[i]][[col]] <- NA
    }
    results[[i]] <- results[[i]][all_cols]
  }
  
  # Add comparison columns
  missing_cols <- setdiff(all_cols, names(comparison_results))
  for(col in missing_cols) {
    comparison_results[[col]] <- NA
  }
  comparison_results <- comparison_results[all_cols]
  
  all_results <- do.call(rbind, results)
  final_results <- cbind(all_results, comparison_results[rep(1, nrow(all_results)), , drop = FALSE])
  
  return(list(results = final_results, models = models))
}

fit_GLM_one = function(af.site, rd.site, sampData, sampDataFull, formulaString, 
                       cmpAll=NULL, dontReport=NULL, ix, site_info=NULL, makePlot=FALSE, 
                       repName=NULL, saveModel=FALSE, useMixedEffects=FALSE, compareModels=FALSE,
                       randomEffectsType="simple") {
  
  if(compareModels && useMixedEffects) {
    model_comparison <- fit_multiple_models(af.site, rd.site, sampData, formulaString, repName, randomEffectsType)
    
    average_freq <- mean(af.site)
    min_freq <- min(af.site)
    max_freq <- max(af.site)
    
    results <- cbind(
      model_comparison$results,
      data.frame(
        average_freq = average_freq,
        min_freq = min_freq,
        max_freq = max_freq,
        stringsAsFactors = FALSE
      )
    )
    
    plot_data = NULL
    if(makePlot && !is.null(sampDataFull) && "gen" %in% names(sampDataFull) && "trt" %in% names(sampDataFull)){
      plot_data <- list(
        raw = data.frame(
          gen = sampDataFull$gen,
          af = af.site,
          trt = sampDataFull$trt,
          rep = if(!is.null(repName)) sampDataFull[[repName]] else NA,
          trt_rep = if(!is.null(repName)) paste(sampDataFull$trt, sampDataFull[[repName]], sep="_") else sampDataFull$trt
        ),
        site_info = site_info,
        formula = formulaString,
        model_type = "model_comparison"
      )
    }
    
    models_to_save = if(saveModel) model_comparison$models else NULL
    
    return(list(results = results, plot_data = plot_data, model = models_to_save))
  } else {
    
    # Original fit_GLM_one logic for single models
    Neff = rd.site
    
    cts = cbind(round(as.numeric(Neff)*as.numeric(af.site)), 
                round(as.numeric(Neff)*(1-as.numeric(af.site))))
    
    if(useMixedEffects) {
      model <- fit_glmer_robust(as.formula(formulaString), sampData, binomial)
      if(is.null(model)) {
        results = data.frame(
          dispersion_parameter = NA,
          average_freq = mean(af.site),
          min_freq = min(af.site),
          max_freq = max(af.site),
          converged = FALSE,
          singular = TRUE,
          overdispersion = NA,
          stringsAsFactors = FALSE
        )
        return(list(results = results, plot_data = NULL, model = NULL))
      }
      
      conv_check <- check_convergence(model)
      overdispersion <- check_overdispersion(model)
      converged <- conv_check$converged
      singular <- conv_check$singular
    } else {
      model = glm(as.formula(formulaString), family="quasibinomial", data=sampData)
      overdispersion <- check_overdispersion(model)
      converged <- TRUE
      singular <- FALSE
    }
    
    residual_deviance <- deviance(model)
    degrees_of_freedom <- df.residual(model)
    dispersion_parameter <- residual_deviance / degrees_of_freedom
    
    average_freq <- mean(af.site)
    min_freq <- min(af.site)
    max_freq <- max(af.site)

    if(useMixedEffects) {
      fixed_effects <- fixef(model)
      fixed_se <- sqrt(diag(vcov(model)))
      fixed_z <- fixed_effects / fixed_se
      fixed_p <- 2 * pnorm(abs(fixed_z), lower.tail = FALSE)
      
      cp <- cbind(fixed_effects, fixed_p)
      rownames(cp) <- names(fixed_effects)

    } else {
      cp = extract_coef_pval(model, cmpAll, dontReport, includeIntercept=TRUE)
    }
    
    results = data.frame(
      matrix(cp[,1], nrow=1, dimnames=list(NULL, paste0("coef.", rownames(cp)))),
      matrix(cp[,2], nrow=1, dimnames=list(NULL, paste0("p.", rownames(cp)))),
      dispersion_parameter = dispersion_parameter,
      average_freq = average_freq,
      min_freq = min_freq,
      max_freq = max_freq,
      converged = converged,
      singular = singular,
      overdispersion = overdispersion,
      stringsAsFactors = FALSE
    )
    

    plot_data = NULL
    if(makePlot && !is.null(sampDataFull) && "gen" %in% names(sampDataFull) && "trt" %in% names(sampDataFull)){
      plot_data <- list(
        raw = data.frame(
          gen = sampDataFull$gen,
          af = af.site,
          trt = sampDataFull$trt,
          rep = if(!is.null(repName)) sampDataFull[[repName]] else NA,
          trt_rep = if(!is.null(repName)) paste(sampDataFull$trt, sampDataFull[[repName]], sep="_") else sampDataFull$trt
        ),
        site_info = site_info,
        formula = formulaString,
        model_type = if(useMixedEffects) "mixed_effects" else "fixed_effects"
      )
    }
    
    model_to_save = if(saveModel) model else NULL
    
    return(list(results = results, plot_data = plot_data, model = model_to_save))
  }
}

process_chunk <- function(chunk_indices, 
                         sampIX, samps, 
                         model.vars, 
                         repName, 
                         repNamecontrast, 
                         cmpAll, 
                         effect1, 
                         effect2, 
                         effect1type, 
                         effect2type, 
                         inter,
                         repNameinter,
                         repNameeff, 
                         dontReport, 
                         args) {

  sampData = set_up_sampData(sampIX, samps, model.vars, repName, repNamecontrast, 
                            cmpAll, effect1, effect1type, effect2, effect2type,
                            trtContrast=args$trtContrast, trtPopGroups=(args$trtPopGroups == "TRUE"))
  
  sampDataFull = samps[sampIX,]
  
  useMixedEffects <- args$mixedEffects == "TRUE"
  compareModels <- args$compareModels == "TRUE"
  randomEffectsType <- args$randomEffectsType
  
  if(args$trtPopGroups == "TRUE" && "trt_pop" %in% names(sampData)){
    if(inter == "TRUE") {
      formulaString = paste0("gen + trt_pop + gen*trt_pop")
    } else {
      formulaString = paste0("gen + trt_pop")
    }
    cat("Using trt_pop model:", formulaString, "\n")
  } else if(useMixedEffects) {
    # IMPROVED: Use simpler random effects structure
    if(randomEffectsType == "simple") {
      if(inter == "TRUE") {
        formulaString = paste0("gen * trt + (1 | ", repName, ")")
      } else {
        formulaString = paste0("gen + trt + (1 | ", repName, ")")
      }
    } else if(randomEffectsType == "nested") {
      if(inter == "TRUE") {
        formulaString = paste0("gen * trt + (1 | ", repName, ":trt)")
      } else {
        formulaString = paste0("gen + trt + (1 | ", repName, ":trt)")
      }
    } else { # crossed
      if(inter == "TRUE") {
        formulaString = paste0("gen * trt + (1 | ", repName, ") + (1 | trt)")
      } else {
        formulaString = paste0("gen + trt + (1 | ", repName, ") + (1 | trt)")
      }
    }
    cat("Using mixed-effects model with", randomEffectsType, "random effects:", formulaString, "\n")
  } else {
    # Original logic for non-trt_pop, non-mixed-effects models
    if(inter == "TRUE" & repNameinter == "TRUE") {
      formulaString = paste0(c(paste0(model.vars, collapse = " + "), 
                              paste(effect1, effect2, sep = "*"), 
                              paste(repName, repNameeff, sep = "*")), collapse = " + ")
    } else if(inter == "TRUE" & repNameinter == "FALSE") { 
      formulaString = paste0(c(paste0(model.vars, collapse = " + "), 
                              paste(effect1, effect2, sep = "*")), collapse = " + ")
    } else if (inter == "FALSE" & repNameinter == "TRUE") {
      formulaString = paste0(c(paste0(model.vars, collapse = " + "), 
                              paste(repName, repNameeff, sep = "*")), collapse = " + ")
    } else if (inter == "FALSE" & repNameinter == "FALSE") {
      formulaString = paste0(model.vars, collapse = " + ")
    }
  }
  
  formulaString = paste0("cts ~ ", formulaString)
  
  makePlot <- args$makePlots == "TRUE"
  saveModel <- args$saveModels == "TRUE"
  
  chunk_formula <- formulaString
  
  results_list <- mclapply(chunk_indices, function(ix) {
    if(ix%%1000 == 0){cat("working on site ",ix,"\n")}
    
    site_info <- data.frame(chrom = sites$chrom[ix], pos = sites$pos[ix])
    
    fit_output <- fit_GLM_one(afmat[ix,sampIX], rd[ix,sampIX], sampData, sampDataFull,
                              formulaString, cmpAll, dontReport, ix, site_info, makePlot, repName, saveModel, useMixedEffects, compareModels, randomEffectsType)
    
    site_key <- paste(sites$chrom[ix], sites$pos[ix], sep="_")
    return(list(result = fit_output$results, plot_data = fit_output$plot_data, model = fit_output$model,
                site_key = site_key))
  }, mc.cores=args$nCores)
  
  results <- do.call(rbind, lapply(results_list, function(x) x$result))
  
  plots <- NULL
  if(makePlot){
    plots <- lapply(results_list, function(x) x$plot_data)
    names(plots) <- lapply(results_list, function(x) x$site_key)
    plots <- plots[!sapply(plots, is.null)]
  }
  
  models <- NULL
  if(saveModel){
    models <- lapply(results_list, function(x) x$model)
    names(models) <- lapply(results_list, function(x) x$site_key)
    models <- models[!sapply(models, is.null)]
  }
  
  print(colnames(results))

  return(list(results = results, plots = plots, models = models, formula = chunk_formula))
}

##########
## MAIN
##########
args <- suppressWarnings(parse_cl_args(commandArgs(trailingOnly=TRUE)))
registerDoMC(cores=args$nCores)
cat("RUNNING calc_GLM.R with the following parameters:\n")
print(args)

# Initialize plot and model lists
all_plots <- list()
all_models <- list()
all_formulas <- list()

# load af data
cat("loading HAFs\n")
load(args$HAFs)  ## should contain sites, samps, afmat
print(samps)

## get chroms and sites
chroms=unique(sites$chrom)
if(!is.na(args$chrom)){
  siteIX=which(sites$chrom==args$chrom)
} else {
  siteIX=1:nrow(sites)
}
if(!is.na(args$testNsites)){
  siteIX=sample(siteIX,args$testNsites)
} 

# set up samp info
samps$cage=samps[[args$repName]]
samps$cage=factor(samps$cage)
cage_set=sort(levels(samps$cage)) 
if(!is.na(args$dropRep)){cage_set=cage_set[cage_set!=args$dropRep]}
sampIX=which(samps$cage %in% cage_set);

# NEW: Filter by Selected Treatments (e.g. "B,T")
if (!is.na(args$selectTrts)) {
  target_trts <- trimws(unlist(strsplit(args$selectTrts, ",")))
  cat("Filtering to keep only treatments:", paste(target_trts, collapse=", "), "\n")
  
  if (!("trt" %in% colnames(samps))) stop("'trt' column not found in samps object.")
  
  trt_matches <- which(samps$trt %in% target_trts)
  sampIX <- intersect(sampIX, trt_matches)
  
  if(length(sampIX) == 0) stop("No samples remained after filtering for treatments: ", args$selectTrts)
}

# set up read depth // effective coverage info
rd=set_up_depthData(args,sites,samps)

# calc GLM
# MODIFIED: Updated variable selection logic to handle trt_pop
if(args$trtPopGroups == "TRUE"){
  # For trt_pop models, we need gen and trt_pop (but trt_pop is created in set_up_sampData)
  variables = c("gen", "trt", args$repName)  # We need these to create trt_pop
  cat("Using trt_pop approach - will create combined treatment-replicate groups\n")
} else {
  # Original logic
  if(args$repNameignore == "TRUE"){
      if(!is.null(args$mainEffect) ){ 
        variables = c(args$mainEffect)
      } else if (!is.na(args$effect2)) {   
        variables = c(args$effect1, args$effect2)
      } else if (!is.na(args$effect1)) { 
        variables = c(args$effect1) 
      }
  } else if (args$repNameignore == "FALSE") {
      if(!is.null(args$mainEffect) ){ 
        variables = c(args$mainEffect, args$repName)
      } else if (!is.na(args$effect2)) {   
        variables = c(args$effect1, args$effect2, args$repName)
      } else if (!is.na(args$effect1)) { 
        variables = c(args$effect1, args$repName) 
      }
  }
}

cat("calculating GLM:",variables,
    "\nusing",length(sampIX),"samples from replicates:",cage_set,"\n")
cat("Random effects type:", args$randomEffectsType, "\n")

filename=paste0(args$outDir,".",args$saveAs)

system.time({
  chunk_size <- 1000  # Smaller chunks for mixed effects models
  num_chunks <- ceiling(length(siteIX) / chunk_size)
  filename <- paste0(args$outDir, ".", args$saveAs)

  for (i in 1:num_chunks) {
    cat("\n=== Processing chunk", i, "of", num_chunks, "===\n")
    chunk_start <- (i - 1) * chunk_size + 1
    chunk_end <- min(i * chunk_size, length(siteIX))
    chunk_indices <- siteIX[chunk_start:chunk_end]
        
    # MODIFIED: Don't filter out coefficients when using trt_pop groups
    dontReport_param <- if(args$trtPopGroups == "TRUE") NULL else args$repName
    
    chunk_output <- process_chunk(
      chunk_indices, 
      sampIX, 
      samps, 
      model.vars = variables, 
      repName = args$repName,
      repNamecontrast = args$repNamecontrast, 
      cmpAll = args$mainEffect, 
      effect1 = args$effect1, 
      effect2 = args$effect2, 
      effect1type = args$effect1type, 
      effect2type = args$effect2type, 
      inter = args$interaction,
      repNameinter = args$repNameinter,
      repNameeff = args$repNameeff, 
      dontReport = dontReport_param, 
      args = args)

    # Extract results, plots, and models
    chunk_results <- chunk_output$results
    
    # Collect plots if they were created
    if(args$makePlots == "TRUE" && !is.null(chunk_output$plots)){
      all_plots <- c(all_plots, chunk_output$plots)
    }
    
    # Collect models if they were saved
    if(args$saveModels == "TRUE" && !is.null(chunk_output$models)){
      all_models <- c(all_models, chunk_output$models)
    }
    
    # Store formula 
    all_formulas[[paste0("chunk", i)]] <- list(
    chunk = i,
    sites = chunk_indices,
    formula = chunk_output$formula,
    chunk_start = chunk_start, 
    chunk_end = chunk_end
    )

    # Combine with chrom and pos columns
    combined_results <- cbind(sites[chunk_indices, c("chrom", "pos")], chunk_results)

    # Write results to file
    if (i == 1) {
      write.table(combined_results, file = filename, sep = ",", row.names = FALSE, col.names = TRUE)
    } else {
      file_conn <- file(filename, open = "a")
      write.table(combined_results, file = file_conn, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
      close(file_conn)
    }
    
    cat("Chunk", i, "completed. Results written to", filename, "\n")
  }
})

# Save environment including plots, models, and formulas
analysis_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
environment_file <- paste0(args$outDir, "_", analysis_timestamp, ".RData")

save.image(file = environment_file)