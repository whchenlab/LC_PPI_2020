## -- data cleansing -- 
require(tidyverse);


## -- modeling, ROC and helper libraries 
require(ROCR);
require(randomForestSRC); ## use randomforestSRC instead ... 
require(caret);
require(pROC);

## -- system helper libs --
require(progress);
require(tictoc); ## show elapsed time

## -- libraries for parallel computing --;
require( doSNOW );
## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
rf.setdata <- function( feat.data, meta.data, grouping = "Group", control = "Healthy" ){

  ## -- elapsed time --
  tic( "Set data " );
  
  ## -- create a list to contain all information ...
  rf <- list();
  
  ## ========================================================================
  ## -- sanity check --
  ## -- 1. if meta.data contains info for all samples;
  sample_names.all <- rownames( feat.data );
  sample_names.not_in_meta_data <- sample_names.all[ !sample_names.all %in% rownames( meta.data ) ];
  if( length( sample_names.not_in_meta_data ) > 0 ){
    stop( "San check failed: the following samples in the feature data are not in the metadata: ",  paste( sample_names.not_in_meta_data, collapse = ", " ), "." );
  }
  
  ## -- 2. if relative abundances --
  feat.range <- range( feat.data, na.rm = TRUE );
  if( feat.range[1] < 0 | feat.range[2] > 1 ){
    warning( "San check warning: feature data may not be relative abundances (ranging from 0 to 1): ", "min = ", feat.range[1], ", max = ", feat.range[2], "\n" );
  }
  
  ## -- 3. if grouping column exists --
  if( ! grouping %in% colnames( meta.data ) ){
    stop( "San check failed: meta data does not contain column with name of '", grouping, "'." );
  }
  
  ## -- 4. if the grouping column contain not two values --
  groups <- table( meta.data[, grouping] );
  control.count <- groups[ control ];
  if( length( table( meta.data[, grouping] ) ) != 2 ){
    stop( "San check failed: the grouping column in meta data should contain EXACTLY TWO categories." );
  }
  
  if( is.na( control.count ) ){
    stop( "San check failed: the control group '", control , "' does not seem to exist in your meta data." );
  }
  
  ## ======================================================================== 
  ## -- set data --
  ## -- 1. back up all raw data --
  rf[["feat.raw"]] <- feat.data;
  rf[["meta.raw"]] <- meta.data;
  
  ## -- 2. add grouping data to feat.data --
  dat <- as.data.frame(feat.data);
  dat[, "Group"] <- as.character( meta.data[ rownames( feat.data ), grouping ] );
  
  ## -- 3. get case and controls 
  group.names <- names(groups);
  case <- group.names[ group.names != control ];
  
  rf[["feat.data"]] <- dat;
  rf[["control.group"]] <- control;
  rf[["case.group"]] <- case;
  rf[["grouping.column"]] <- grouping;
  
  ## ======================================================================== 
  ## -- wrap up and return rf --
  cat( "All look great.", "\n");
  
  cat( "-----------------------------------------------------------", "\n");
  cat( "  In total, the feature data contains ", nrow(feat.data), " samples, with ", ncol(feat.data), " features.\n", sep = "" );
  cat( "  Case group: '", case, "', w/ ", groups[case], " samples, \n", sep = "");
  cat( "  Control group: '", control, "', w/ ", groups[control], " samples, \n", sep = "");
  cat( "-----------------------------------------------------------", "\n\n");
  
  cat( " Data size: ", format( object.size(rf), units = "MB" ), "\n", sep = "");
  
  ## -- 
  toc(); ## show elapsed time 
  ## -- return --
  return( rf );
}
## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------

rf.train <- function( rf, fold = 10, resample = 10, parallel = F, num.cores = 10, class.weights = NULL ){
  
  ## -- requires ... 
  require(purrr); 
  
  ## -- elapsed time --
  tic( "Model training " );
  
  ## ========================================================================
  ## -- sanity check --
  ## 1. check if input contains necessary data 
  if( sum( c("feat.data", "control.group", "case.group") %in% names(rf)  ) < 3  ){
    stop( "San check failed: input data does not contain all required slots: feat.data, control.group, case.group, please try from the begining." );
  }
  
  
  ## ========================================================================
  ## -- get dataa from the input list --
  dat <- rf[["feat.data"]]; 
  dat$Group <- as.factor( dat$Group );
  
  meta.data <- rf[["meta.raw"]];
  grouping <- rf[["grouping.column"]]; ## the column in meta.data for the grouping information, i.e. CASE, CONTROL
  control.group <- rf[["control.group"]];
  case.group <- rf[["case.group"]];
  
  ## -- 
  cat("Training models on ", nrow( dat ), " samples, with ", ncol( dat ) - 1, " features.\n", sep = "" );
  cat("Case group: ", case.group, ", control group: ", control.group, "\n", sep = "" );
  cat("  - running ", resample, " resample and ", fold , " fold modeling ... \n", sep = "");
  
  ## -- 
  ## -- class weights --
  if( isTRUE( class.weights ) ){ 
      cat("  - using class weight option of randomForest classifier ... \n" );
  } else {
      cat("  - NOT using class weight option of randomForest classifier ... \n" );
  }
  
  ## ========================================================================
  ## -- modeling --
  
  ##require group name is "Group"
  ## -- create confusion matrix --
  sample.split <- createMultiFolds( as.vector( dat$Group ), times = resample, k = fold  );
  
  ## -- create a few lists to save final data 
  rf.models <- list();
  pred.matrix <- list();
  confusion.result <- list();
  
  ## --- ----
  if( parallel ){
    ## ========================================================================
    ##    parallel version 
    ## ========================================================================
    ## 1. decide number of cpu cores to use 
    available.cores <- parallel::detectCores() - 1;
    if( num.cores > available.cores  ){
      num.cores = available.cores;
    }
    
    if( num.cores < 1 ){
      num.cores = 1;
    }
    
    cat( "  - running parallel modeling using ", num.cores, " cpu cores ... \n", sep = "" );
    
    ## -- 2. modeling using designated CPU cores  --
    cl <- makeSOCKcluster(num.cores);
    registerDoSNOW(cl);
    
    pb <- txtProgressBar(min=1, max=length( sample.split ), style=3);
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    
    model.results <- 
      foreach( split = sample.split, .packages="randomForestSRC", .options.snow=opts ) %dopar% {
        ## -- get data ready --
        train.data <- dat[ split , ];
        test.data <- dat[ -split , ];
        
        ## -- train, and get model prediction results ... 
        if(  isTRUE( class.weights ) ){
            this.model <- rfsrc( Group ~ .,train.data, proximity = T, importance = T, ntree = 500,
                                 case.wt = randomForestSRC:::make.wt( train.data$Group ),
                                 sampsize = randomForestSRC:::make.size( train.data$Group ));
        } else {
            this.model <- rfsrc(Group ~ .,train.data, proximity = T, importance = T, ntree = 500);
        }
        ## 
        pred <- predict( this.model, test.data, type = "prob");
        this.pred.matrix <- pred$predicted ;
        rownames( this.pred.matrix ) <- rownames( test.data );
        
        ## 
        this.confusion.result <- predict( this.model, test.data, type= "response");
        
        return( list( "model" = this.model,
                      "pred.matrix" = this.pred.matrix,
                      "confusion.result" = this.confusion.result) );
      }  ## end of foreach ... 
    close(pb);
    stopCluster(cl);
    
    ## -- 3. post process to get required data ... 
    iters <- names( sample.split );
    for( idx in 1:length(model.results) ){
      iter <- iters[ idx ];
      result.list <- model.results[[idx]];
      
      ## -- get data from list and save them to other lists ... 
      rf.models[[iter]] <- result.list[["model"]];
      pred.matrix[[iter]] <- result.list[["pred.matrix"]];
      confusion.result[[iter]] <- result.list[["confusion.result"]];
    }
    
  } else {
    ## ========================================================================
    ##    regular version 
    ## ========================================================================
    ## -- make a progress bar for modeling -- 
    pb1 <- progress_bar$new( format = "Modeling [:bar] :percent in :elapsed", clear = FALSE, total =  fold * resample );
    
    ## -- do modeling ... 
    for (iter in names(sample.split)) {
      ## -- get data ready --
      train.data <- dat[sample.split[[iter]], ]
      test.data <- dat[-sample.split[[iter]], ]
      
      ## -- train, and predict ... 
      if(  isTRUE( class.weights ) ){
          this.model <- rfsrc( Group ~ .,train.data, proximity = T, importance = T, ntree = 500,
                               case.wt = randomForestSRC:::make.wt(train.data$Group),
                               sampsize = randomForestSRC:::make.size(train.data$Group));
      } else {
          this.model <- rfsrc(Group ~ .,train.data, proximity = T, importance = T, ntree = 500);
      }
      
      ## -- predict --
      pred <- predict( this.model, test.data, type = "prob");
      this.pred.matrix <- pred$predicted;
      rownames( this.pred.matrix ) <- rownames( test.data );
      
      ## -- confusion matrix --
      this.confusion.result <- predict(this.model, test.data, type= "response")
      
      ## -- save results ...
      rf.models[[iter]] = this.model;
      pred.matrix[[iter]] <- this.pred.matrix;
      confusion.result[[iter]] <- this.confusion.result;
      
      ## -- tick --
      pb1$tick();
    }
  }
  
  ## --- importance scores --- ##
  importance <- lapply(rf.models, function(data){ imp <-  data$importance[, c("all")] * 100 }); ## -- ... 
  importance.df <- purrr::reduce(importance, cbind)
  importance.df <- data.frame(importance.df, stringsAsFactors = F);
  names(importance.df) <- names(rf.models);
  
  
  ## -- assemble rf --
  ## -- modeling results ...
  rf[["sample.split"]] <- sample.split;
  rf[["rf.models"]] <- rf.models;
  rf[["pred.matrix"]] <- pred.matrix;
  rf[["confusion.result"]] <- confusion.result;
  rf[["importance.df"]] <- importance.df;
  
  cat( " Data size: ", format( object.size(rf), units = "MB" ), "\n", sep = "");
  
  ## -- return ... 
  toc(); ## show elapsed time 
  ## -- return --
  return( rf );
}
## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
rf.evaluate <- function( rf ){
  
  ## -- requires ... 
  require(tictoc); ## show elapsed time
  require(tidyverse);
  
  ## -- modeling, ROC and helper libraries 
  require(ROCR);
  require(pROC);
  
  ## -- elapsed time --
  tic( "Model performance evaluation " );
  
  ## ========================================================================
  ## -- sanity check --
  ## 1. check if input contains necessary data 
  if( sum( "pred.matrix" %in% names(rf) ) < 1  ){
    stop( "San check failed: input data does not contain required slot: pred.matrix, please check if you have trained the models." );
  }
  
  if( sum( "meta.raw" %in% names(rf) ) < 1  ){
    stop( "San check failed: input data does not contain required slot: meta.raw, please check if you have the meta data correctly assigned." );
  }
  
  if( sum( "grouping.column" %in% names(rf) ) < 1  ){
    stop( "San check failed: input data does not contain required slot: grouping.column, please check if you have the meta data correctly assigned." );
  }
  
  ## -- 
  cat("Evaluating model performance ... \n", sep = "");
  
  ## ========================================================================
  ## -- get data from previous results --
  pred.matrix       <- rf[["pred.matrix"]];
  meta.data         <- rf[["meta.raw"]];
  grouping          <- rf[["grouping.column"]]; ## the column in meta.data for the grouping information, i.e. CASE, CONTROL
  control.group     <- rf[["control.group"]];
  case.group        <- rf[["case.group"]];
  
  ## ========================================================================
  ## -- get model performance ... 
  ## -- 1. iterate the prediction matrix, combine all predictions, and calculate mean scores for each sample --
  res.pred.df <- pred.matrix %>% 
    map_df( function(x) { tibble( Sample = rownames(x), Control = x[, control.group], Case = x[, case.group ] ) } ) %>% 
    group_by( Sample ) %>% 
    summarize( Control = mean(Control), Case = mean(Case) ) %>% as.data.frame( . );
  
  colnames( res.pred.df ) <- c("Sample", control.group, case.group );
  
  ## -- 2. summarize prediction results and add the original grouping data ... 
  res.pred.df$Predicted <- apply(res.pred.df, 1, function(data){ names(which.max(data[ c( control.group, case.group ) ]))} );
  res.pred.df$Group <-  meta.data[ res.pred.df$Sample, grouping ];
  
  ## -- 3. calculate AUC & CI --
  rf.new <- evaluate.internal.use( res.pred.df, control.group, case.group );
  
  toc(  ); ## show elapsed time ... 
  return( rf.new );
}

## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
rf.external.validation <- function( rf, ext.feat.data, ext.meta.data ){
    
    ## -- elapsed time --
    tic( "External validation " );
    
    ## ========================================================================
    ## -- sanity check & run predictions--
    ## 1. check if input contains necessary data 
    if( sum( c("rf.models", "control.group", "case.group", "grouping.column" ) %in% names( rf )  ) < 4  ){
        stop( "San check failed: input data does not contain all required slots: 'rf.models', 'control.group', 'case.group', 'grouping.column', please check your input." );
    }
    
    ## -- get data from tf ... 
    grouping <- rf[["grouping.column"]]; ## the column in meta.data for the grouping information, i.e. CASE, CONTROL
    control.group <- rf[["control.group"]];
    case.group <- rf[["case.group"]];
    
    ## -- 2. try to validate feat and meta data by build a new rf object using rf.setdata --
    rf.ext.data <- rf.setdata( ext.feat.data, ext.meta.data, grouping = grouping, control = control.group );
    
    ## -- 3. get feat data with and without groups, and use rf.predict to do predictions --
    feat.w.grp <- rf.ext.data$feat.data;
    feat.w.o.grp <- feat.w.grp[ , !colnames(feat.w.grp) %in% "Group" ];
    
    ext.data.pred <- rf.predict( rf, feat.w.o.grp );
    
    ## -- 4. evaluate performance --
    res.pred.df <- ext.data.pred$predicted.results.df;
    res.pred.df$Group <-  ext.meta.data[ res.pred.df$Sample, grouping ];
    
    #### ======= THE FOLLOWING CODES ARE THE SAME FROM rf.evaluate =======
    rf.new <- evaluate.internal.use( res.pred.df, control.group, case.group );
    
    toc(  ); ## show elapsed time ... 
    return( rf.new );
}

## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
evaluate.internal.use <- function( res.pred.df, control.group, case.group ){
    
    forestpred <- ROCR::prediction( res.pred.df[, case.group ], res.pred.df[, "Group" ],
                                    label.ordering = c( control.group, case.group ));
    
    forestperf <- ROCR::performance(forestpred, "tpr", "fpr");
    data.roc.df <- data.frame( fpr = forestperf@x.values[[1]],  tpr = forestperf@y.values[[1]] );
    
    ## -- 4. use pROC to calculate AUC with CI -
    auc <- pROC::roc( res.pred.df[, "Group"], res.pred.df[, case.group ], levels = c( control.group, case.group ), ci = TRUE, plot = FALSE );
    auc.ci.values <- round( as.numeric( auc$ci ), digits = 2 ); ## a vector of three
    auc.string <- paste0( auc.ci.values[2], " (95%CI: ",  auc.ci.values[1], "-",  auc.ci.values[3], ")")
    
    ## -- 5. plot --
    data.roc.plot <- 
        ggplot(data.roc.df, aes(x = fpr, y = tpr) ) + geom_line() + theme_bw() + 
        labs(x = "False Positve Rate", y = "True Positive Rate", title = "ROC plot")  + 
        annotate( "text", x = .75, y = .75, label = paste0( "AUC:", auc.string ), size = 4) +
        theme( plot.title = element_text(hjust = .5, face = "bold", size = 11), 
               axis.text = element_text(size = 9) );
    
    print( data.roc.plot );
    
    ## ========================================================================
    ### calculate confusion matrix --
    confusion <- table( res.pred.df[,  c( "Group", "Predicted" ) ] );
    errors <- round( 1 - diag( confusion ) / rowSums( confusion ), digits = 3 );
    errors.overall <- round(  1 - ( sum(diag(confusion)) / sum( confusion ) ) , digits = 3 );
    
    ## ========================================================================
    ## -- print on-screen statistics ...  
    cat("\n--------------------------------------------\n");
    cat( "Overall performance (AUC): ", auc.string, "\n\n" );
    cat( "Predicted vs. Truth (Group): \n" );
    cat( "Control group: ", control.group, "\n" );
    cat( "Case group: ", case.group, "\n\n" );
    print( confusion );
    cat( "---------\n" );
    cat( "Error rates: \n" );
    print( errors );
    cat("\nOverall error rate: ", errors.overall, "\n", sep = "");
    cat("--------------------------------------------\n\n");
    ## ========================================================================
    ## -- 
    
    rf.new <- list( "roc.df" = data.roc.df, 
                    "roc.plot" = data.roc.plot, 
                    "auc" = auc.ci.values[2], 
                    "auc.with.ci" = auc.ci.values, 
                    "auc.obj" = auc, 
                    "auc.string" = auc.string,
                    "pred.df" = res.pred.df, 
                    "error.rates" = errors,
                    "overall.error.rate" = errors.overall);
    
    cat( "Data size: ", format( object.size(rf.new), units = "MB" ), "\n", sep = "");
    cat("---\n");
    
    return(rf.new);
}
## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------

