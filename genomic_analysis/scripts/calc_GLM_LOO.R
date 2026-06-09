###############
## Libraries
###############
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(doMC)))
suppressWarnings(suppressMessages(library(argparser)))
suppressWarnings(suppressMessages(library(lme4)))

# Control settings
glmerControl_settings <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))

## FUNCTIONS
#################
parse_cl_args=function(args){
  p <- arg_parser("Calculate GLM coefficients with Leave-One-Out Analysis")
  p <- add_argument(p, "HAFs", help=".Rdata file")
  p <- add_argument(p, "--readDepth", help=".RDS file", default=NA)
  p <- add_argument(p, "--effectiveCov", help="cov", default=NA)
  p <- add_argument(p, "--chrom", help="chrom", default=NA)
  p <- add_argument(p, "--poolSize", help="N", default=100,type="integer")
  p <- add_argument(p, "--mainEffect", help="var", default="NA")
  p <- add_argument(p, "--effect1", help="var", default="NA")
  p <- add_argument(p, "--effect1type", help="type", default="NA")  
  p <- add_argument(p, "--effect2", help="var", default="NA")
  p <- add_argument(p, "--effect2type", help="type", default="NA")
  p <- add_argument(p, "--interaction", help="bool", default="FALSE")
  p <- add_argument(p, "--repName", help="col", default="rep")
  p <- add_argument(p, "--repNamecontrast", help="contrast", default="NA")
  p <- add_argument(p, "--repNameignore", help="bool", default="FALSE")
  p <- add_argument(p, "--repNameinter", help="bool", default="FALSE")
  p <- add_argument(p, "--repNameeff", help="var", default="NA")
  p <- add_argument(p, "--testNsites", help="N", default=NA,type="integer")
  p <- add_argument(p, "--nCores", help="N", default=2,type="integer")
  p <- add_argument(p, "--saveAs", help="fmt", default="csv")
  p <- add_argument(p, "--outDir", help="dir", default=".")
  p <- add_argument(p, "--trtContrast", help="ref", default=NA)
  p <- add_argument(p, "--selectTrts", help="trts", default=NA)
  p <- add_argument(p, "--trtPopGroups", help="bool", default="FALSE")
  p <- add_argument(p, "--mixedEffects", help="bool", default="FALSE")
  p <- add_argument(p, "--randomEffectsType", help="type", default="simple")

  args <- parse_args(p)
  
  if (args$mainEffect == "NA") {args$mainEffect=NULL}
  if (args$effect1 == "NA") {args$effect1=NA}
  if (args$effect2 == "NA") {args$effect2=NA}
  
  return(args)
}

set_up_sampData=function(sampIX, samps, model.vars, repName, repNamecontrast, cmpAll, effect1, effect1type, effect2, effect2type, trtContrast, trtPopGroups=FALSE){
  sampData=samps[sampIX,]%>%dplyr::select(all_of(model.vars))
  sampData[] <- lapply(sampData, function(x) if(is.factor(x)) droplevels(x) else x)

  if(trtPopGroups && "trt" %in% names(sampData) && repName %in% names(sampData)){
    sampData$trt_pop <- interaction(sampData$trt, sampData[[repName]], sep="_", drop=TRUE)
  }
  
  if(!is.na(trtContrast) && "trt" %in% names(sampData)){
    sampData$trt = factor(sampData$trt)
    all_levels = levels(sampData$trt)
    new_levels = c(trtContrast, all_levels[all_levels != trtContrast])
    sampData$trt = factor(sampData$trt, levels = new_levels)
  }
    
  if(!is.na(effect1)){
    if(effect1type == "continuous") sampData[, effect1] = as.numeric(sampData[, effect1])
    else if(effect1type == "categorical") sampData[, effect1] = as.factor(sampData[, effect1])
  }
    
  if(!is.na(effect2)){
    if(effect2type == "continuous") sampData[, effect2] = as.numeric(sampData[, effect2])
    else if(effect2type == "categorical") sampData[, effect2] = as.factor(sampData[, effect2])
  }

  if(!is.na(repName) && repNamecontrast == "sum"){
      sampData[, repName] = as.factor(sampData[, repName])
      contrasts(sampData[, repName]) <- contr.sum(length(levels(sampData[, repName])))
  }
 
  return(sampData)
} 

set_up_depthData=function(args,sites,samps){
  if(!is.na(args$readDepth)){
    rd=data.frame(readRDS(args$readDepth))
    ss=match(samps$SourceFile,colnames(rd))
    rd=rd[,ss]
  } else if(!is.numeric(args$effectiveCov)){
    df.ec=fread(args$effectiveCov);
    df.ec$chrom <- factor(df.ec$chrom)
    mat.ec = t(xtabs(ec ~ sampID + chrom, df.ec))
    cc=match(sites$chrom,rownames(mat.ec))
    ss=match(samps$SourceFile,colnames(mat.ec))
    rd=round(mat.ec[cc,ss])
  } else{
    rd=matrix(args$effectiveCov,nrow(sites),nrow(samps))
  }
  return(rd)
}


fit_model_internal <- function(formulaString, data, family="quasibinomial", useMixed=FALSE, returnModel=FALSE) {
    model <- NULL
    try({
        if(useMixed){
            model <- glmer(as.formula(formulaString), data=data, family=binomial, control=glmerControl_settings)
        } else {
            model <- glm(as.formula(formulaString), family=family, data=data)
        }
    }, silent=TRUE)
    
    if(returnModel) return(model)
    return(NULL)
}

fit_GLM_LOO = function(af.site, rd.site, sampData, formulaString, 
                       repName, useMixedEffects, trtPopGroups, 
                       inter_term_name="gen:trt") {
  
  Neff = rd.site
  cts = cbind(round(as.numeric(Neff)*as.numeric(af.site)), 
              round(as.numeric(Neff)*(1-as.numeric(af.site))))

  modelData <- cbind(data.frame(cts=I(cts)), sampData)
  

  full_model <- fit_model_internal(formulaString, modelData, useMixed=useMixedEffects, returnModel=TRUE)
  

  full_res <- list(coef=NA, pval=NA, converged=FALSE)
  target_coef_name <- NULL
  
  if(!is.null(full_model)) {
    full_res$converged <- if(useMixedEffects) is.null(full_model@optinfo$conv$lme4$messages) else TRUE
    
    if(useMixedEffects) {
        coefs <- fixef(full_model)
        se <- sqrt(diag(vcov(full_model)))
        z <- coefs / se
        p <- 2 * pnorm(abs(z), lower.tail = FALSE)
        res_table <- cbind(coefs, p)
    } else {
        res_table <- summary(full_model)$coefficients[,c(1,4),drop=FALSE]
    }
    
    term_idx <- grep(":", rownames(res_table))
    if(length(term_idx) > 0) {
        target_coef_name <- rownames(res_table)[tail(term_idx, 1)]
    } else {
        target_coef_name <- rownames(res_table)[nrow(res_table)]
    }
    
    full_res$coef <- res_table[target_coef_name, 1]
    full_res$pval <- res_table[target_coef_name, 2]
  }


  trts <- unique(as.character(modelData$trt))
  
  loo_definitions <- list()
  

  if(length(trts) == 2) {
      reps1 <- sort(unique(as.character(modelData[[repName]][modelData$trt == trts[1]])))
      reps2 <- sort(unique(as.character(modelData[[repName]][modelData$trt == trts[2]])))
      

      combo_grid <- expand.grid(r1 = reps1, r2 = reps2, stringsAsFactors = FALSE)
      

      for(r in 1:nrow(combo_grid)) {
          drop_list <- list()
          drop_list[[trts[1]]] <- combo_grid$r1[r]
          drop_list[[trts[2]]] <- combo_grid$r2[r]
          loo_definitions[[length(loo_definitions) + 1]] <- drop_list
      }
      
  } else {
      reps_by_trt <- list()
      for(t in trts) {
          reps_by_trt[[t]] <- sort(unique(as.character(modelData[[repName]][modelData$trt == t])))
      }
      n_loo <- min(sapply(reps_by_trt, length))
      for(i in 1:n_loo) {
          drop_list <- lapply(reps_by_trt, function(x) x[i])
          loo_definitions[[length(loo_definitions) + 1]] <- drop_list
      }
  }
  
  loo_stats <- list()
  
  for(i in 1:length(loo_definitions)) {
      reps_to_drop <- loo_definitions[[i]]

      drop_mask <- rep(FALSE, nrow(modelData))
      for(t in names(reps_to_drop)) {
          drop_rep <- reps_to_drop[[t]]
          drop_mask <- drop_mask | (modelData$trt == t & modelData[[repName]] == drop_rep)
      }
      
      train_data <- modelData[!drop_mask, ]
      test_data <- modelData[drop_mask, ]

      loo_model <- fit_model_internal(formulaString, train_data, useMixed=useMixedEffects, returnModel=TRUE)
      
      loo_coef <- NA
      loo_rmse <- NA
      loo_pred_mean <- NA
      
      if(!is.null(loo_model)) {
          if(useMixedEffects) {
              c_val <- fixef(loo_model)
              if(!is.null(target_coef_name) && target_coef_name %in% names(c_val)) {
                  loo_coef <- c_val[target_coef_name]
              } else {
                   loo_coef <- c_val[length(c_val)] 
              }
          } else {
              c_val <- coef(loo_model)
              if(!is.null(target_coef_name) && target_coef_name %in% names(c_val)) {
                  loo_coef <- c_val[target_coef_name]
              } else {
                   loo_coef <- c_val[length(c_val)]
              }
          }
          
          try({
              preds <- predict(loo_model, newdata=test_data, type="response", allow.new.levels=TRUE)
              actuals <- test_data$cts[,1] / (test_data$cts[,1] + test_data$cts[,2])
              loo_rmse <- sqrt(mean((preds - actuals)^2, na.rm=TRUE))
              loo_pred_mean <- mean(preds, na.rm=TRUE)
          }, silent=TRUE)
      }
      
      loo_stats[[paste0("pair_", i)]] <- list(
          coef = loo_coef,
          rmse = loo_rmse,
          pred = loo_pred_mean
      )
  }
  
  vec_coef <- sapply(loo_stats, function(x) x$coef)
  vec_rmse <- sapply(loo_stats, function(x) x$rmse)
  
  results <- data.frame(
      full_coef = full_res$coef,
      full_pval = full_res$pval,
      full_converged = full_res$converged,

      loo_mean_coef = mean(vec_coef, na.rm=TRUE),
      loo_sd_coef = sd(vec_coef, na.rm=TRUE),
      loo_cv_coef = sd(vec_coef, na.rm=TRUE) / abs(mean(vec_coef, na.rm=TRUE)),
      loo_mean_rmse = mean(vec_rmse, na.rm=TRUE),
      
      loo_n_success = sum(!is.na(vec_coef)),
      
      stringsAsFactors = FALSE
  )
  
  for(i in 1:length(loo_stats)) {
      results[[paste0("loo",i,"_coef")]] <- loo_stats[[i]]$coef
      results[[paste0("loo",i,"_rmse")]] <- loo_stats[[i]]$rmse
  }
  
  return(results)
}

process_chunk_loo <- function(chunk_indices, sampIX, samps, variables, repName, args) {
  
  sampData = set_up_sampData(sampIX, samps, variables, repName, args$repNamecontrast, 
                            args$mainEffect, args$effect1, args$effect1type, args$effect2, args$effect2type,
                            trtContrast=args$trtContrast, trtPopGroups=(args$trtPopGroups == "TRUE"))
  
  useMixedEffects <- args$mixedEffects == "TRUE"
  
  inter <- args$interaction
  if(args$trtPopGroups == "TRUE"){
      formulaString = if(inter=="TRUE") "gen * trt_pop" else "gen + trt_pop"
  } else if(useMixedEffects) {
      randomPart <- if(args$randomEffectsType == "nested") paste0("(1 | ", repName, ":trt)") else paste0("(1 | ", repName, ")")
      formulaString = if(inter=="TRUE") paste0("gen * trt + ", randomPart) else paste0("gen + trt + ", randomPart)
  } else {
      if(inter == "TRUE") formulaString = paste0(paste(variables, collapse=" + "), " + ", args$effect1, "*", args$effect2)
      else formulaString = paste0(variables, collapse=" + ")
  }
  
  formulaString = paste0("cts ~ ", formulaString)

  results_list <- mclapply(chunk_indices, function(ix) {
    
    fit_res <- fit_GLM_LOO(afmat[ix,sampIX], rd[ix,sampIX], sampData, 
                           formulaString, repName, useMixedEffects, args$trtPopGroups)
    
    return(fit_res)
  }, mc.cores=args$nCores)
  
  results <- do.call(rbind, results_list)
  return(results)
}

##########
## MAIN
##########
args <- suppressWarnings(parse_cl_args(commandArgs(trailingOnly=TRUE)))
registerDoMC(cores=args$nCores)
cat("RUNNING calc_GLM_LOO.R with Combinatorial Pairing\n")

load(args$HAFs) # loads sites, samps, afmat

siteIX <- if(!is.na(args$chrom)) which(sites$chrom==args$chrom) else 1:nrow(sites)
if(!is.na(args$testNsites)) siteIX <- sample(siteIX, args$testNsites)

samps$cage = samps[[args$repName]]
cage_set = sort(unique(samps$cage))
sampIX = which(samps$cage %in% cage_set)

if (!is.na(args$selectTrts)) {
  target_trts <- trimws(unlist(strsplit(args$selectTrts, ",")))
  sampIX <- intersect(sampIX, which(samps$trt %in% target_trts))
}

rd = set_up_depthData(args, sites, samps)

if(args$trtPopGroups == "TRUE"){
  variables = c("gen", "trt", args$repName) 
} else {
  variables = c(args$mainEffect, args$repName) # Simplified defaults
  if(!is.na(args$effect1)) variables <- c(variables, args$effect1)
  if(!is.na(args$effect2)) variables <- c(variables, args$effect2)
  variables <- unique(variables[!is.na(variables)])
}

chunk_size <- 500 # Smaller chunks due to high compute load
num_chunks <- ceiling(length(siteIX) / chunk_size)
filename <- paste0(args$outDir, ".", args$saveAs)

for (i in 1:num_chunks) {
  chunk_start <- (i - 1) * chunk_size + 1
  chunk_end <- min(i * chunk_size, length(siteIX))
  chunk_indices <- siteIX[chunk_start:chunk_end]
  
  chunk_res <- process_chunk_loo(chunk_indices, sampIX, samps, variables, args$repName, args)
  
  combined_results <- cbind(sites[chunk_indices, c("chrom", "pos")], chunk_res)

  if (i == 1) {
    write.table(combined_results, file = filename, sep = ",", row.names = FALSE, col.names = TRUE)
  } else {
    file_conn <- file(filename, open = "a")
    write.table(combined_results, file = file_conn, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
    close(file_conn)
  }
}
