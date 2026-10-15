# =============================================================================
# EMF_statistics.py 的纯 R 转换，并重建「视觉修订版」的五幅图
# 使用：在 RStudio 中以 UTF-8 打开，点击 Source。
# 不调用 Python。随机森林使用 ranger；R/NumPy/ sklearn 的随机实现不同，
# 因此Bootstrap区间、RF数值和空间置换p值不保证逐位等于Python。
# 原Python脚本没有绘图函数，本脚本的图形依据同目录PNG与配色规范重建。
# =============================================================================
INPUT_FILE <- Sys.getenv('EMF_INPUT_FILE',
  unset='D:/Users/Lenovo/Desktop/原始/3.2/3.2.2/样点+EMF.csv')
OUTPUT_DIR <- Sys.getenv('EMF_OUTPUT_DIR',
  unset='D:/Users/Lenovo/Desktop/原始/3.2/3.2.3/R_视觉复现版')
SEED <- 20261014L
BOOTSTRAP_B <- 2000L
OUTER_K <- 5L
INNER_K <- 4L
PERMUTATION_REPEATS <- 100L
TUNING_TREES <- 300L
FINAL_TREES <- 1000L
RF_THREADS <- 2L
DPI <- 600L
AUTO_INSTALL <- TRUE
# 可选：填写本脚本生成的EMF_complete_results.rds，仅重画图，默认完整重算。
REUSE_RESULTS <- Sys.getenv('EMF_REUSE_STATS')

main <- function() {
  options(stringsAsFactors=FALSE,timeout=300)
  mm <- paste(R.version$major,strsplit(R.version$minor,'.',fixed=TRUE)[[1]][1],sep='.')
  custom_lib <- Sys.getenv('EMF_R_LIBRARY')
  persistent_lib <- file.path(dirname(OUTPUT_DIR),'EMF_R_library')
  user_lib <- file.path(Sys.getenv('LOCALAPPDATA'),'R','win-library',mm)
  paths <- c(custom_lib,persistent_lib,user_lib,.libPaths())
  .libPaths(paths[nzchar(paths)&dir.exists(paths)])
  pkgs <- c('ranger','officer','rvg','jsonlite')
  missing <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
  if (length(missing)) {
    if (!AUTO_INSTALL) stop('请安装：',paste(missing,collapse=', '))
    lib <- if(nzchar(custom_lib)) custom_lib else persistent_lib
    dir.create(lib,recursive=TRUE,showWarnings=FALSE); .libPaths(c(lib,.libPaths()))
    install.packages(missing,repos='https://cloud.r-project.org',lib=lib,
                    type=if(.Platform$OS.type=='windows') 'binary' else 'source')
  }
  if (!'min.bucket' %in% names(formals(ranger::ranger))) stop('请更新ranger到支持min.bucket的版本。')
  out <- OUTPUT_DIR
  if(file.exists(file.path(out,'EMF_statistics.txt'))&&!nzchar(REUSE_RESULTS))
    out <- paste0(out,'_',format(Sys.time(),'%Y%m%d_%H%M%S'))
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  csvdir <- file.path(out,'统计明细CSV'); dir.create(csvdir,showWarnings=FALSE)
  vars <- c('Carbon','Nutrient','Enzyme','Vegetation'); av <- c(vars,'EMF')
  original <- c('Plot ID','Region','Longitude','Latitude','EMF-carbon',
                'EMF-nutrient','EMF-enzyme','EMF-vegetation','EMF-avg')
  data <- read.csv(INPUT_FILE,check.names=FALSE,fileEncoding='UTF-8-BOM')
  norm <- function(x) gsub('[[:space:]_-]','',tolower(trimws(x)))
  if(ncol(data)!=9 || !identical(norm(names(data)),norm(original))) stop('请检查A–I列名称及顺序。')
  names(data) <- c('PlotID','Region','Longitude','Latitude',av)
  for(j in 3:9) data[[j]] <- suppressWarnings(as.numeric(data[[j]]))
  if(anyNA(data)||anyDuplicated(data$PlotID)||!all(is.finite(as.matrix(data[,3:9]))))
    stop('原Python脚本要求无缺失、无重复Plot ID及无非有限数值。请先核查数据。')
  d <- data; d$Region <- factor(d$Region)
  X <- as.matrix(d[,vars]); y <- d$EMF; n <- length(y)
  regions <- split(seq_len(n),d$Region)
  if(min(lengths(regions))<max(OUTER_K,INNER_K)) stop('区域样本数不足当前折数。')
  mean_error <- max(abs(y-rowMeans(X)))
  if(mean_error>=1e-8) stop('EMF-avg不是四维等权均值，不能直接应用本次指数方差分解。')
  md5 <- unname(tools::md5sum(INPUT_FILE))
  tables <- list()
  save_table <- function(name,tab) {
    tables[[name]] <<- tab
    write.csv(tab,file.path(csvdir,paste0(name,'.csv')),row.names=FALSE,fileEncoding='UTF-8')
  }
  stars <- function(q) ifelse(q<.001,'***',ifelse(q<.01,'**',ifelse(q<.05,'*','')))
  form <- reformulate(vars,'EMF')
  metrics <- function(obs,pr) c(R2=1-sum((obs-pr)^2)/sum((obs-mean(obs))^2),
    RMSE=sqrt(mean((obs-pr)^2)),MAE=mean(abs(obs-pr)))
  settings <- list(seed=SEED,bootstrap_B=BOOTSTRAP_B,outer_K=OUTER_K,inner_K=INNER_K,
    permutations=PERMUTATION_REPEATS,tuning_trees=TUNING_TREES,final_trees=FINAL_TREES,
    outer_repeats=1L,importance='Unweighted mean of fold-specific held-out Delta R2',
    min_node_size=1L,bootstrap='Within Region; fixed region sample counts')
  if(nzchar(REUSE_RESULTS)) {
    res <- readRDS(REUSE_RESULTS)
    if(!identical(md5,res$source_md5)||!identical(settings,res$settings))
      stop('复用结果的数据或统计设置不匹配。')
    tables <- res$tables
    for(nm in names(tables)) save_table(nm,tables[[nm]])
    message('已核验数据与设置，仅重画图。')
  } else {
    message('相关、回归和指数构造核验')
    save_table('01_data_audit',data.frame(check=c('n','regions','missing','duplicate_ID',
      'equal_mean_max_error','source','MD5'),value=c(n,paste(names(regions),lengths(regions),collapse=', '),
      0,0,format(mean_error,digits=16),INPUT_FILE,md5)))
    save_table('02_descriptive',do.call(rbind,lapply(av,function(v) data.frame(variable=v,n=n,
      mean=mean(d[[v]]),SD=sd(d[[v]]),min=min(d[[v]]),Q25=unname(quantile(d[[v]],.25)),
      median=median(d[[v]]),Q75=unname(quantile(d[[v]],.75)),max=max(d[[v]])))))
    save_table('03_region_descriptive',do.call(rbind,lapply(names(regions),function(g)
      do.call(rbind,lapply(av,function(v) { z<-d[[v]][regions[[g]]]
        data.frame(Region=g,variable=v,n=length(z),mean=mean(z),SD=sd(z),min=min(z),max=max(z))})))))
    corrs <- list(); target <- list()
    for(method in c('pearson','spearman')) {
      pairs <- do.call(rbind,lapply(combn(av,2,simplify=FALSE),function(a) {
        args <- list(x=d[[a[1]]],y=d[[a[2]]],method=method)
        if(method=='spearman') args$exact <- FALSE
        ct <- do.call(cor.test,args)
        data.frame(method=method,var1=a[1],var2=a[2],n=n,r=unname(ct$estimate),p=ct$p.value,
          CI95_low=if(method=='pearson') ct$conf.int[1] else NA_real_,
          CI95_high=if(method=='pearson') ct$conf.int[2] else NA_real_)
      }))
      pairs$q_BH_10pairs <- p.adjust(pairs$p,'BH')
      rr <- diag(1,5); pp <- qq <- matrix(NA_real_,5,5)
      dimnames(rr)<-dimnames(pp)<-dimnames(qq)<-list(av,av)
      for(i in seq_len(nrow(pairs))) {
        a<-pairs$var1[i]; b<-pairs$var2[i]
        rr[a,b]<-rr[b,a]<-pairs$r[i]; pp[a,b]<-pp[b,a]<-pairs$p[i]
        qq[a,b]<-qq[b,a]<-pairs$q_BH_10pairs[i]
      }
      save_table(paste0('04_',method,'_r_matrix'),data.frame(variable=av,rr,row.names=NULL))
      save_table(paste0('05_',method,'_p_matrix'),data.frame(variable=av,pp,row.names=NULL))
      save_table(paste0('06_',method,'_q_matrix'),data.frame(variable=av,qq,row.names=NULL))
      save_table(paste0('07_',method,'_pairs'),pairs)
      zz<-pairs[pairs$var2=='EMF',];zz$q_BH_4tests<-p.adjust(zz$p,'BH');target[[method]]<-zz
      corrs[[method]]<-list(r=rr,p=pp,q=qq,pairs=pairs)
    }
    save_table('08_dimension_EMF_correlations',do.call(rbind,target))
    sensitivity <- do.call(rbind,lapply(c('All',names(regions),'Adjusted for Region'),function(g) {
      ids <- if(g %in% c('All','Adjusted for Region')) seq_len(n) else regions[[g]]
      do.call(rbind,lapply(vars,function(v) {
        a<-d[[v]][ids]; b<-y[ids]; other<-rowMeans(X[ids,setdiff(vars,v),drop=FALSE]); df<-length(ids)-2
        if(g=='Adjusted for Region') {
          a<-resid(lm(a~d$Region)); b<-resid(lm(b~d$Region)); other<-resid(lm(other~d$Region))
          df<-n-nlevels(d$Region)-1
        }
        do.call(rbind,lapply(c('EMF','Other three'),function(tg) {
          r<-cor(a,if(tg=='EMF') b else other); t<-r*sqrt(df/(1-r*r))
          data.frame(analysis=g,dimension=v,target=tg,n=length(ids),df=df,r=r,p=2*pt(abs(t),df,lower.tail=FALSE))
        }))
      }))
    }))
    sensitivity$q_BH_4tests <- ave(sensitivity$p,interaction(sensitivity$analysis,sensitivity$target),FUN=function(p)p.adjust(p,'BH'))
    save_table('09_region_and_other_three_correlations',sensitivity)
    fit<-lm(form,d); sm<-summary(fit); cf<-sm$coefficients; ci<-confint(fit)
    A<-model.matrix(fit); e<-resid(fit); h<-hatvalues(fit); bread<-solve(crossprod(A))
    vcHC3<-bread%*%crossprod(A,A*as.numeric((e/(1-h))^2))%*%bread
    seHC3<-sqrt(diag(vcHC3)); nu<-df.residual(fit); crit<-qt(.975,nu)
    vif<-sapply(vars,function(v)1/(1-summary(lm(reformulate(setdiff(vars,v),v),d))$r.squared))
    coef_table<-data.frame(term=rownames(cf),b=cf[,1],SE=cf[,2],t=cf[,3],p=cf[,4],
      CI95_low=ci[,1],CI95_high=ci[,2],HC3_SE=seHC3,
      HC3_p=2*pt(abs(coef(fit)/seHC3),nu,lower.tail=FALSE),
      HC3_CI95_low=coef(fit)-crit*seHC3,HC3_CI95_high=coef(fit)+crit*seHC3,
      beta_standardized=c(NA,coef(fit)[vars]*apply(X,2,sd)/sd(y)),VIF=c(NA,vif),row.names=NULL)
    save_table('10_OLS_coefficients_and_HC3',coef_table)
    bp<-n*summary(lm(e^2~X))$r.squared
    save_table('11_OLS_fit_and_diagnostics',data.frame(n=n,R2=sm$r.squared,adjusted_R2=sm$adj.r.squared,
      F=unname(sm$fstatistic[1]),F_p=pf(sm$fstatistic[1],4,nu,lower.tail=FALSE),
      RMSE=sqrt(mean(e^2)),MAE=mean(abs(e)),AIC=AIC(fit),BIC=BIC(fit),
      df_model=4,df_residual=nu,condition_number_raw=kappa(A,exact=TRUE),
      condition_number_standardized=kappa(cbind(1,scale(X)),exact=TRUE),
      Shapiro_p=if(n<=5000)shapiro.test(e)$p.value else NA_real_,Breusch_Pagan_p=pchisq(bp,4,lower.tail=FALSE),
      note='Rounding-level residuals; these diagnostics and p values are not ecological inference.'))
    save_table('12_OLS_plot_diagnostics',data.frame(PlotID=d$PlotID,Region=d$Region,observed=y,
      fitted=fitted(fit),residual=e,leverage=h,Cook_distance=cooks.distance(fit)))
    capture.output(sm,coef_table,file=file.path(out,'OLS_raw_output.txt'))
    # 同Python：加权中心协方差的代数计算，yy默认为重构的四维度均值。
    partition<-function(xx,yy=rowMeans(xx),weights=rep(1/nrow(xx),nrow(xx))) {
      weights<-weights/sum(weights); xc<-sweep(xx,2,colSums(xx*weights),'-')
      yc<-yy-sum(yy*weights); cc<-crossprod(xc,xc*weights); cy<-colSums(xc*as.numeric(weights*yc))
      vy<-sum(weights*yc^2); shares<-.25*cy/vy; own<-.0625*diag(cc)/vy
      r2<-numeric(16)
      for(mask in 1:15) {
        idx<-which(as.logical(intToBits(mask)[1:4])); C<-cc[idx,idx,drop=FALSE]
        eg<-eigen(C,symmetric=TRUE); tol<-max(eg$values)*1e-15
        inv<-eg$vectors%*%diag(ifelse(eg$values>tol,1/eg$values,0),nrow=length(idx))%*%t(eg$vectors)
        r2[mask+1]<-as.numeric(t(cy[idx])%*%inv%*%cy[idx]/vy)
      }
      shap<-sapply(1:4,function(j) {
        total<-0
        for(mask in 0:15) { bits<-as.logical(intToBits(mask)[1:4])
          if(!bits[j]) { k<-sum(bits); w<-factorial(k)*factorial(3-k)/factorial(4)
            total<-total+w*(r2[mask+2^(j-1)+1]-r2[mask+1]) }
        };total
      })
      list(variance=shares,Shapley=shap,own=own,subsets=r2)
    }
    allocation<-partition(X); stopifnot(abs(sum(allocation$variance)-1)<1e-9,abs(sum(allocation$Shapley)-1)<1e-9)
    message('区域分层Bootstrap：',BOOTSTRAP_B,'次')
    set.seed(SEED); bs<-bl<-matrix(NA_real_,BOOTSTRAP_B,4)
    for(b in seq_len(BOOTSTRAP_B)) {
      ids<-unlist(lapply(regions,function(ii)sample(ii,length(ii),replace=TRUE)))
      z<-partition(X[ids,,drop=FALSE]);bs[b,]<-z$variance;bl[b,]<-z$Shapley
    }
    colnames(bs)<-colnames(bl)<-vars; bci<-apply(bs,2,quantile,c(.025,.975));lci<-apply(bl,2,quantile,c(.025,.975))
    contribution<-data.frame(dimension=vars,index_weight=.25,covariance_share=allocation$variance,
      own_variance_share=allocation$own,allocated_covariance_share=allocation$variance-allocation$own,
      covariance_CI95_low=bci[1,],covariance_CI95_high=bci[2,],Shapley_R2=allocation$Shapley,
      Shapley_share=allocation$Shapley/sum(allocation$Shapley),Shapley_CI95_low=lci[1,],Shapley_CI95_high=lci[2,],
      bootstrap_top_frequency_covariance=tabulate(max.col(bs,ties.method='first'),4)/BOOTSTRAP_B,
      bootstrap_top_frequency_Shapley=tabulate(max.col(bl,ties.method='first'),4)/BOOTSTRAP_B,row.names=NULL)
    save_table('13_contribution_partitioning',contribution)
    save_table('14_all_bootstrap_draws',data.frame(draw=1:BOOTSTRAP_B,
      setNames(as.data.frame(bs),paste0('CovShare_',vars)),setNames(as.data.frame(bl),paste0('Shapley_',vars))))
    save_table('15_all_subset_R2',data.frame(mask=0:15,predictors=sapply(0:15,function(m) {
      ids<-which(as.logical(intToBits(m)[1:4]));if(length(ids))paste(vars[ids],collapse='+') else '(Intercept only)'
    }),R2=allocation$subsets))
    csens<-do.call(rbind,lapply(c(names(regions),'Region balanced'),function(g) {
      ids<-if(g=='Region balanced')seq_len(n) else regions[[g]]
      weights<-if(g=='Region balanced')1/as.numeric(table(d$Region)[d$Region]) else rep(1,length(ids))
      z<-partition(X[ids,,drop=FALSE],weights=weights)
      data.frame(analysis=g,dimension=vars,n=length(ids),covariance_share=z$variance,Shapley_R2=z$Shapley)
    }))
    save_table('16_region_and_balanced_contributions',csens)
    diff<-do.call(rbind,lapply(combn(4,2,simplify=FALSE),function(ij) do.call(rbind,lapply(c('covariance','Shapley'),function(method) {
      B<-if(method=='covariance')bs else bl; estimate<-if(method=='covariance')allocation$variance else allocation$Shapley
      delta<-B[,ij[1]]-B[,ij[2]]
      data.frame(method=method,dimension1=vars[ij[1]],dimension2=vars[ij[2]],difference=estimate[ij[1]]-estimate[ij[2]],
        CI95_low=unname(quantile(delta,.025)),CI95_high=unname(quantile(delta,.975)),bootstrap_frequency_greater=mean(delta>0))
    }))))
    save_table('17_bootstrap_contribution_differences',diff)

    # 对应sklearn: max_features .5/1 -> mtry 2/4；min_samples_leaf -> min.bucket。
    # ranger停止分裂的阈值与sklearn含义不同。min.node.size=1允许2个样本的节点分裂，
    # 对应sklearn默认min_samples_split=2；叶片下限由min.bucket单独约束。
    # 不调用Python，也不把旧Python结果当作本次R结果。
    folds<-function(group,k,seed) {
      set.seed(seed);f<-integer(length(group))
      for(ii in split(seq_along(group),group))f[sample(ii)]<-rep(seq_len(k),length.out=length(ii));f
    }
    grid<-expand.grid(mtry=c(2L,4L),min.bucket=c(1L,3L,5L),max.depth=c(0L,6L))
    grow<-function(rows,par,trees,seed,importance='impurity')ranger::ranger(
      x=d[rows,vars,drop=FALSE],y=y[rows],num.trees=trees,mtry=par$mtry,min.node.size=1L,
      min.bucket=par$min.bucket,max.depth=par$max.depth,replace=TRUE,sample.fraction=1,
      splitrule='variance',importance=importance,num.threads=RF_THREADS,write.forest=TRUE,seed=seed)
    pred<-function(model,nd)as.numeric(predict(model,data=nd,num.threads=RF_THREADS)$predictions)
    tune<-function(rows,k,seed,forest_seed) {
      ff<-folds(d$Region[rows],k,seed); candidates<-grid
      candidates$mean_test_MSE<-candidates$sd_test_MSE<-NA_real_
      for(j in seq_len(nrow(grid))) {
        scores<-sapply(seq_len(k),function(v) {
          tr<-rows[ff!=v];te<-rows[ff==v]
          f<-grow(tr,grid[j,],TUNING_TREES,forest_seed,'none')
          mean((y[te]-pred(f,d[te,vars,drop=FALSE]))^2)
        })
        # sklearn GridSearchCV：各折分数等权均值、总体标准差。
        candidates$mean_test_MSE[j]<-mean(scores)
        candidates$sd_test_MSE[j]<-sqrt(mean((scores-mean(scores))^2))
      }
      candidates$rank<-rank(candidates$mean_test_MSE,ties.method='min')
      list(best=candidates[which.min(candidates$mean_test_MSE),1:3,drop=FALSE],grid=candidates)
    }
    message('RF：一次5折外层/4折内层嵌套CV，每维度每外层折100次置换')
    outer<-folds(d$Region,OUTER_K,SEED);oof<-linear_oof<-numeric(n)
    perf<-raw_perm<-fold_imp<-tuning<-list()
    for(k in seq_len(OUTER_K)) {
      tr<-which(outer!=k);te<-which(outer==k)
      chosen<-tune(tr,INNER_K,SEED+k,SEED+k)
      tuning[[k]]<-cbind(outer_fold=k,chosen$grid)
      f<-grow(tr,chosen$best,FINAL_TREES,SEED+k)
      nd<-d[te,vars,drop=FALSE];pr<-pred(f,nd);oof[te]<-pr
      linear_oof[te]<-predict(lm(form,d[tr,]),newdata=d[te,])
      ss<-metrics(y[te],pr)
      perf[[k]]<-data.frame(fold=k,n_train=length(tr),n_test=length(te),R2=ss[1],RMSE=ss[2],MAE=ss[3],
        inner_best_RMSE=sqrt(min(chosen$grid$mean_test_MSE)),chosen$best,row.names=NULL)
      set.seed(SEED+1000L+k)
      for(j in seq_along(vars)) {
        # 一次批量预测全部置换，减少反复启动线程的开销。
        xp<-nd[rep(seq_len(nrow(nd)),PERMUTATION_REPEATS),,drop=FALSE]
        for(b in seq_len(PERMUTATION_REPEATS))xp[(b-1)*nrow(nd)+seq_len(nrow(nd)),j]<-sample(nd[[j]])
        pp<-matrix(pred(f,xp),nrow=length(te),ncol=PERMUTATION_REPEATS)
        den<-sum((y[te]-mean(y[te]))^2)
        delta<-(colSums((pp-y[te])^2)-sum((pr-y[te])^2))/den
        raw_perm[[length(raw_perm)+1]]<-data.frame(fold=k,dimension=vars[j],repeat_id=1:PERMUTATION_REPEATS,delta_R2=delta)
        mdi<-f$variable.importance/sum(f$variable.importance)
        fold_imp[[length(fold_imp)+1]]<-data.frame(fold=k,dimension=vars[j],
          permutation_delta_R2_mean=mean(delta),permutation_repeat_SD=sd(delta),MDI=mdi[j])
      }
      message('外层折 ',k,'/',OUTER_K,'，R²=',sprintf('%.4f',ss[1]))
    }
    fs<-tune(seq_len(n),OUTER_K,SEED+100L,SEED)
    full<-grow(seq_len(n),fs$best,FINAL_TREES,SEED)
    fi<-do.call(rbind,fold_imp);mdi<-full$variable.importance/sum(full$variable.importance)
    top_per_fold<-sapply(1:OUTER_K,function(k){z<-fi[fi$fold==k,];z$dimension[which.max(z$permutation_delta_R2_mean)]})
    imp<-do.call(rbind,lapply(seq_along(vars),function(j) {
      z<-fi[fi$dimension==vars[j],]
      data.frame(dimension=vars[j],permutation_delta_R2=mean(z$permutation_delta_R2_mean),
        fold_SD=sd(z$permutation_delta_R2_mean),MDI_full=mdi[j],MDI_fold_mean=mean(z$MDI),
        MDI_fold_SD=sd(z$MDI),top_fold_count=sum(top_per_fold==vars[j]),row.names=NULL)
    }))
    imp$permutation_relative_positive_share<-pmax(imp$permutation_delta_R2,0)/sum(pmax(imp$permutation_delta_R2,0))
    save_table('18_RF_importance_summary',imp);save_table('19_RF_outer_fold_performance',do.call(rbind,perf))
    save_table('20_RF_outer_fold_importance',fi);save_table('21_RF_all_permutations',do.call(rbind,raw_perm))
    save_table('22_RF_nested_tuning',do.call(rbind,tuning));save_table('23_RF_final_tuning',fs$grid)
    score<-metrics(y,oof);lm_score<-metrics(y,linear_oof)
    rfmetrics<-data.frame(n=n,nested_OOF_R2=score[1],nested_OOF_RMSE=score[2],nested_OOF_MAE=score[3],
      train_R2=metrics(y,pred(full,d[,vars]))[1],OOB_R2=metrics(y,full$predictions)[1],
      linear_OOF_R2=lm_score[1],linear_OOF_RMSE=lm_score[2],outer_folds=OUTER_K,inner_folds=INNER_K,
      grid_combinations=12,permutations=PERMUTATION_REPEATS,row.names=NULL)
    save_table('24_RF_model_performance',rfmetrics)
    parameters<-c(list(backend='ranger',num.trees=FINAL_TREES,min.node.size=1L,replace=TRUE,sample.fraction=1,
      splitrule='variance',importance='normalized impurity',seed=SEED),as.list(fs$best[1,]))
    save_table('25_RF_parameters',data.frame(parameter=names(parameters),value=vapply(parameters,as.character,character(1))))
    oofrows<-data.frame(PlotID=d$PlotID,Region=d$Region,Longitude=d$Longitude,Latitude=d$Latitude,
      outer_fold=outer,observed_EMF=y,RF_OOF_prediction=oof,RF_OOF_residual=y-oof,OLS_OOF_prediction=linear_oof)
    save_table('26_plot_OOF_predictions',oofrows)
    transfer<-transfer_imp<-list()
    for(g in names(regions)) {
      tr<-setdiff(seq_len(n),regions[[g]]);te<-regions[[g]]
      f<-grow(tr,data.frame(mtry=4L,min.bucket=3L,max.depth=0L),FINAL_TREES,SEED)
      nd<-d[te,vars,drop=FALSE];pr<-pred(f,nd);ss<-metrics(y[te],pr)
      transfer[[g]]<-data.frame(train_region=paste(setdiff(names(regions),g),collapse='+'),test_region=g,
        n_train=length(tr),n_test=length(te),R2=ss[1],RMSE=ss[2],MAE=ss[3],row.names=NULL)
      set.seed(SEED+2000L)
      for(j in seq_along(vars)) {
        xp<-nd[rep(seq_len(nrow(nd)),PERMUTATION_REPEATS),,drop=FALSE]
        for(b in 1:PERMUTATION_REPEATS)xp[(b-1)*nrow(nd)+seq_len(nrow(nd)),j]<-sample(nd[[j]])
        pp<-matrix(pred(f,xp),nrow=length(te));den<-sum((y[te]-mean(y[te]))^2)
        delta<-(colSums((pp-y[te])^2)-sum((pr-y[te])^2))/den
        transfer_imp[[length(transfer_imp)+1]]<-data.frame(test_region=g,dimension=vars[j],delta_R2=mean(delta),permutation_SD=sd(delta))
      }
    }
    save_table('27_region_transfer_performance',do.call(rbind,transfer))
    save_table('28_region_transfer_importance',do.call(rbind,transfer_imp))
    # 与原Python相同的全局8NN逆距离邻域，以及全局残差随机置换。
    rad<-pi/180;lon<-d$Longitude*rad;lat<-d$Latitude*rad
    aa<-sin(outer(lat,lat,'-')/2)^2+outer(cos(lat),cos(lat))*sin(outer(lon,lon,'-')/2)^2
    distance<-6371*2*asin(sqrt(pmin(pmax(aa,0),1)));diag(distance)<-Inf;W<-matrix(0,n,n)
    for(i in seq_len(n)){ix<-order(distance[i,])[1:8];W[i,ix]<-1/pmax(distance[i,ix],1e-6)}
    W<-W/rowSums(W);W<-(W+t(W))/2
    moran<-function(v){z<-v-mean(v);n/sum(W)*as.numeric(t(z)%*%W%*%z)/sum(z^2)}
    obs<-moran(y-oof);set.seed(SEED+3000L);null<-replicate(999,moran(sample(y-oof)))
    spatial<-data.frame(variable='RF OOF residual',Moran_I=obs,permutations=999,
      p_upper_tail=(1+sum(null>=obs))/1000,seed=SEED+3000L,
      note='Global 8NN inverse haversine distance; exploratory; no spatial correction of other tests.')
    save_table('29_spatial_residual_diagnostic',spatial)
    res<-list(source_md5=md5,settings=settings,tables=tables,corr=corrs,contribution=contribution,
      bootstrap_variance=bs,bootstrap_Shapley=bl,RF_importance=imp,RF_metrics=rfmetrics,
      predictions=oofrows,transfer=do.call(rbind,transfer),RF_model=full,linear_model=fit,data=d)
  }
  saveRDS(res,file.path(out,'EMF_complete_results.rds'),compress='xz')
  jsonlite::write_json(list(source=INPUT_FILE,MD5=md5,settings=settings,tables=tables),
    file.path(out,'EMF_results.json'),pretty=TRUE,auto_unbox=TRUE,na='null',digits=15)
  con<-file(file.path(out,'EMF_statistics.txt'),'wt',encoding='UTF-8')
  writeLines(c('EMF Python-to-R 统计转换','对应原Python：一次5×4嵌套CV、每折100次置换、2000次分层Bootstrap。',
    '数学构造、方差分解、Shapley点估计与原Python相同。随机估计不要求逐位相同。',
    'EMF-avg为等权四维均值，回归p值与统计贡献不代表独立因果效应。'),con)
  for(nm in names(tables)) if(!nm %in% c('14_all_bootstrap_draws','21_RF_all_permutations','26_plot_OOF_predictions'))
    writeLines(c('',paste0('===== ',nm,' ====='),capture.output(print(tables[[nm]],row.names=FALSE,digits=8))),con)
  close(con);capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
  # 五幅图的原生R/grid绘图。蓝/橙与点形双重编码，数值由本次R结果生成。
  INK<-'#22272B';MUTED<-'#697780';BLUE<-'#0072B2';ORANGE<-'#E69F00'
  GREY<-'#B0B9C1';RULE<-'#C5CED5';GRID<-'#E9EDF0';FONT<-'Arial'
  labels<-c(Carbon='Carbon',Nutrient='Nutrient',Enzyme='Enzyme',Vegetation='Vegetation',EMF='EMF average')
  captions<-c(
    paste0('Figure 1. Pearson correlations among four functional dimensions and their equal-weight EMF average (n = ',n,
      '). Off-diagonal cells show r and Benjamini-Hochberg significance across 10 pairwise tests. The identity diagonal is neutral. ',
      '* q < 0.05, ** q < 0.01, *** q < 0.001. Component-average associations include mathematical coupling. Tests assume independent plots.'),
    paste0('Figure 2. Plot-level relationships of (a) carbon, (b) nutrient, (c) enzyme and (d) vegetation function with EMF average. ',
      'Circles indicate GQ and diamonds indicate FT. Lines show pooled univariate least-squares fits. Annotations report Pearson r ',
      'and unadjusted two-sided p (n = ',n,'). The response contains each dimension, so these associations do not identify independent ecological effects.'),
    paste0('Figure 3. Random forest importance in descending order. (a) Validation R² decrease after ',PERMUTATION_REPEATS,
      ' permutations per feature in each outer fold of Region-stratified nested cross-validation. Thin error bars show the SD of the ',
      OUTER_K,' fold means, not a confidence interval. (b) Normalized impurity importance from the final forest. ',
      'Blue highlights the leading predictor. Importance describes reconstruction of the composite response and depends on predictor covariance.'),
    paste0('Figure 4. Paired estimates of contribution. Orange circles partition index variance as 0.25 Cov(Xj, M) / Var(M), ',
      'where M is the reconstructed EMF average. Blue squares show Shapley (LMG) allocation of regression R² across all 24 predictor entry orders. ',
      'Error bars are percentile 95% intervals from ',BOOTSTRAP_B,' resamples within Region. Each allocation sums to 100%. ',
      'These conditional intervals assume independent plots within Region. The allocations describe statistical variation, not causal contributions.'),
    paste0('Figure 5. (a) Held-out EMF predictions from ',OUTER_K,' outer folds, with nested Region-stratified tuning. ',
      'The dashed line is perfect agreement. Circles indicate GQ and diamonds indicate FT. (b) Region transfer using a separate fixed forest ',
      'with ',FINAL_TREES,' trees, all four predictors per split and minimum leaf size 3. Random folds can retain spatial dependence. ',
      'The response is a deterministic composite of predictors, and transfer results apply only to the sampled Regions.'))
  W<-518.4 # 虚拟排版单位为pt，对应7.2英寸宽；同一布局同时导出TIFF和矢量PPTX。
  heights<-c(450,450,390,400,420)
  canvas<-function(H) {
    grid::grid.newpage()
    S<-grid::convertWidth(grid::unit(1,'npc'),'inches',valueOnly=TRUE)*72/W
    text<-function(s,x,y,size=9.2,col=INK,just='left',bold=FALSE,rot=0) {
      grid::grid.text(s,x=x/W,y=1-y/H,just=c(just,'centre'),rot=rot,
        gp=grid::gpar(fontfamily=FONT,fontsize=size*S,col=col,fontface=if(bold)'bold' else 'plain'))
    }
    line<-function(x1,y1,x2,y2,col=RULE,lwd=.5,lty=1) {
      grid::grid.lines(x=c(x1,x2)/W,y=1-c(y1,y2)/H,
        gp=grid::gpar(col=col,lwd=lwd*S/.75,lty=lty))
    }
    rect<-function(x,y,w,h,fill,border=NA,lwd=.3) {
      grid::grid.rect(x=x/W,y=1-y/H,width=w/W,height=h/H,just=c('left','top'),
        gp=grid::gpar(fill=fill,col=border,lwd=lwd*S/.75))
    }
    point<-function(x,y,col,shape='circle',radius=1.8,alpha=1,border=NA) {
      fill<-grDevices::adjustcolor(col,alpha.f=alpha)
      if(shape=='circle')grid::grid.circle(x=x/W,y=1-y/H,r=grid::unit(radius*S,'pt'),
        gp=grid::gpar(fill=fill,col=border,lwd=.25*S/.75))
      if(shape=='square')rect(x-radius,y-radius,2*radius,2*radius,fill,border)
      if(shape=='diamond')grid::grid.polygon(x=c(x-radius,x,x+radius,x)/W,
        y=1-c(y,y-radius,y,y+radius)/H,gp=grid::gpar(fill=fill,col=border,lwd=.25*S/.75))
    }
    measure<-function(s,size=7.2) {
      grid::pushViewport(grid::viewport(gp=grid::gpar(fontfamily=FONT,fontsize=size*S)))
      val<-grid::convertWidth(grid::stringWidth(s),'inches',valueOnly=TRUE)*72/S
      grid::popViewport();val
    }
    wrap<-function(s,width,size=7.2) {
      words<-strsplit(s,' +')[[1]];lines<-character();current<-''
      for(word in words) {
        candidate<-if(nzchar(current))paste(current,word) else word
        if(measure(candidate,size)>width && nzchar(current)){lines<-c(lines,current);current<-word}
        else current<-candidate
      };c(lines,current)
    }
    caption<-function(s,y) {
      lines<-wrap(s,W-82,7.2)
      for(j in seq_along(lines))text(lines[j],49,y+(j-1)*9.5,7.2,'#526075')
    }
    list(text=text,line=line,rect=rect,point=point,caption=caption)
  }
  header<-function(g,title,subtitle) {
    g$text(title,49,27,13.4);g$text(subtitle,49,44,8.3,MUTED)
  }
  axis<-function(g,left,top,width,height,xlim,ylim,xticks,yticks,
                 xlab='',ylab='',x_digits=0,y_digits=1,grid_x=FALSE) {
    xx<-function(v)left+(v-xlim[1])/diff(xlim)*width
    yy<-function(v)top+height-(v-ylim[1])/diff(ylim)*height
    if(grid_x)for(v in xticks)g$line(xx(v),top,xx(v),top+height,GRID,.4)
    g$line(left,top+height,left+width,top+height,RULE,.6)
    if(length(yticks))g$line(left,top,left,top+height,RULE,.6)
    for(v in xticks) {
      g$line(xx(v),top+height,xx(v),top+height+2.5,MUTED,.4)
      g$text(sprintf(paste0('%.',x_digits,'f'),v),xx(v),top+height+9,9,MUTED,'centre')
    }
    for(v in yticks) {
      g$line(left-2.5,yy(v),left,yy(v),MUTED,.4)
      g$text(sprintf(paste0('%.',y_digits,'f'),v),left-6,yy(v),9,INK,'right')
    }
    if(nzchar(xlab))g$text(xlab,left+width/2,top+height+26,9.5,just='centre')
    if(nzchar(ylab))g$text(ylab,left-36,top+height/2,9.5,just='centre',rot=90)
    list(x=xx,y=yy)
  }
  plot_points<-function(g,x,y,reg,trans) {
    # 先GQ后FT，与参考图的层次保持一致；FT用钻石形且略大。
    for(group in c('GQ','FT',setdiff(unique(as.character(reg)),c('GQ','FT')))) {
      ids<-which(as.character(reg)==group)
      for(i in ids)g$point(trans$x(x[i]),trans$y(y[i]),if(group=='FT')ORANGE else BLUE,
        shape=if(group=='FT')'diamond' else 'circle',radius=if(group=='FT')2.8 else 1.65,
        alpha=if(group=='FT').8 else .6,border=if(group=='FT')'#B3822C' else '#FFFFFF80')
    }
  }
  region_key<-function(g,x,y) {
    g$point(x,y,BLUE,radius=1.7,alpha=.65);g$text('GQ',x+11,y,9)
    g$point(x+45,y,ORANGE,'diamond',2.7,.8,'#B3822C');g$text('FT',x+56,y,9)
  }
  draw1<-function(with_caption=TRUE) {
    H<-if(with_caption)heights[1] else heights[1]-50;g<-canvas(H)
    header(g,'Functional associations',paste0('Pearson correlation   /   ',n,' plots   /   significance adjusted across 10 pairs'))
    r<-res$corr$pearson$r;q<-res$corr$pearson$q
    pal<-grDevices::colorRampPalette(c('#C88D23','#E9D5B0','#FAF9F6','#AECBDD','#155F88'))(257)
    colour<-function(v)pal[1+round((max(-1,min(1,v))+1)/2*256)]
    left<-140;top<-48;cell<-58
    for(i in 1:5) {
      g$text(labels[av[i]],left-8,top+(i-.5)*cell,9.3,just='right',bold=i==5)
      g$text(labels[av[i]],left+(i-.5)*cell,top+5*cell+24,9.1,MUTED,'centre',i==5,rot=35)
      for(j in 1:i) {
        cx<-left+(j-.5)*cell;cy<-top+(i-.5)*cell
        g$rect(left+(j-1)*cell+.4,top+(i-1)*cell+.4,cell-.8,cell-.8,
          if(i==j)'#FAF9F6' else colour(r[i,j]))
        col<-if(i==j)MUTED else if(abs(r[i,j])>.65)'white' else INK
        g$text(sprintf('%.2f',r[i,j]),cx,cy-if(i==j)0 else 3.5,11.2,col,'centre')
        if(i!=j)g$text(stars(q[i,j]),cx,cy+11,7.5,col,'centre')
      }
    }
    g$line(left,top+4*cell,left+5*cell,top+4*cell,'#AAB5BE',.5)
    lx<-342;ly<-104;lw<-124
    g$text('Pearson r',lx,86,8.3,MUTED)
    for(k in 1:128)g$rect(lx+(k-1)*lw/128,ly,lw/128+.03,5.6,pal[1+2*(k-1)])
    for(k in c(-1,0,1)){xx<-lx+(k+1)*lw/2;g$line(xx,ly+5.8,xx,ly+8,MUTED,.35);g$text(k,xx,ly+14,8,MUTED,'centre')}
    g$text('* q < 0.05   ** q < 0.01',lx,147,8,MUTED)
    g$text('*** q < 0.001',lx,158,8,MUTED)
    if(with_caption)g$caption(captions[1],410)
  }
  draw2<-function(with_caption=TRUE) {
    H<-if(with_caption)heights[2] else heights[2]-50;g<-canvas(H)
    header(g,'Dimensions and average multifunctionality',paste0('Pooled least-squares associations   /   FT: ',sum(d$Region=='FT'),
      ' plots   /   GQ: ',sum(d$Region=='GQ'),' plots'))
    region_key(g,418,37)
    target<-res$corr$pearson$pairs;target<-target[target$var2=='EMF',]
    for(j in 1:4) {
      left<-if(j%%2==1)64 else 299;top<-if(j<=2)85 else 237
      g$text(letters[j],left-2,top-17,10,bold=TRUE);g$text(vars[j],left+14,top-17,10)
      pp<-target$p[match(vars[j],target$var1)]
      g$text(paste0('r = ',sprintf('%.3f',target$r[match(vars[j],target$var1)]),
        '   p ',if(pp<.001)'< 0.001' else paste0('= ',sprintf('%.3f',pp))),left+181,top-17,7.8,MUTED,'right')
      xr<-range(X[,j]);xr<-xr+c(-1,1)*diff(xr)*.06
      xt<-if(j==4)seq(-1.5,1.5,.5) else seq(-2,2,1)
      xt<-xt[xt>=xr[1]&xt<=xr[2]]
      tr<-axis(g,left,top,181,98,xr,c(-1.65,1.85),xt,c(-1.5,-.5,.5,1.5),
        paste0(vars[j],' function score'),if(j%%2==1)'EMF average' else '',x_digits=if(j==4)1 else 0)
      f<-lm(y~X[,j]);b<-coef(f);xx<-seq(min(X[,j]),max(X[,j]),length.out=200);yy<-b[1]+b[2]*xx
      ok<-yy>=-1.65&yy<=1.85
      g$line(tr$x(min(xx[ok])),tr$y(yy[which(ok)[1]]),tr$x(max(xx[ok])),tr$y(tail(yy[ok],1)),INK,.9)
      plot_points(g,X[,j],y,d$Region,tr)
    }
    if(with_caption)g$caption(captions[2],409)
  }
  draw3<-function(with_caption=TRUE) {
    H<-if(with_caption)heights[3] else heights[3]-52;g<-canvas(H)
    met<-res$RF_metrics;imp<-res$RF_importance
    header(g,'Random forest importance',sprintf('Nested out-of-fold R² = %.3f   /   RMSE = %.3f   /   %s trees',
      met$nested_OOF_R2,met$nested_OOF_RMSE,format(FINAL_TREES,big.mark=',')))
    for(panel in 1:2) {
      left<-if(panel==1)77 else 311;wid<-171
      g$text(if(panel==1)'a' else 'b',left-13,69,10,bold=TRUE)
      g$text(if(panel==1)'Held-out permutation' else 'Forest impurity',left,69,10)
      g$text(if(panel==1)'Mean ± SD across five folds' else 'Normalized MDI',left,88,8.5,MUTED)
      z<-imp[order(if(panel==1)-imp$permutation_delta_R2 else -imp$MDI_full),]
      value<-if(panel==1)z$permutation_delta_R2 else 100*z$MDI_full
      err<-if(panel==1)z$fold_SD else rep(0,4)
      xmax<-if(panel==1)max(.5,ceiling(max(value+err)*10)/10+.1) else max(50,ceiling((max(value)+3)/10)*10)
      x<-function(v)left+v/xmax*wid
      for(t in if(panel==1)seq(0,xmax,.1) else seq(0,xmax,10)) {
        g$line(x(t),100,x(t),287,GRID,.4);g$text(if(panel==1)sprintf('%.1f',t) else t,x(t),297,9,MUTED,'centre')
      }
      g$line(left,287,left+wid,287,RULE,.6)
      for(i in 1:4) {
        yy<-119+(i-1)*49
        g$text(z$dimension[i],left-4,yy,9,just='right')
        g$rect(left,yy-10.5,x(value[i])-left,21,if(i==1)BLUE else GREY)
        if(panel==1) {
          g$line(x(value[i]-err[i]),yy,x(value[i]+err[i]),yy,INK,.6)
          g$line(x(value[i]-err[i]),yy-2.5,x(value[i]-err[i]),yy+2.5,INK,.6)
          g$line(x(value[i]+err[i]),yy-2.5,x(value[i]+err[i]),yy+2.5,INK,.6)
        }
        g$text(sprintf(if(panel==1)'%.4f' else '%.2f',value[i]),x(value[i]+err[i])+3,yy,8.5)
      }
      g$text(if(panel==1)'Decrease in validation R²' else 'Mean decrease in impurity (%)',left+wid/2,318,9.5,just='centre')
    }
    if(with_caption)g$caption(captions[3],348)
  }
  draw4<-function(with_caption=TRUE) {
    H<-if(with_caption)heights[4] else heights[4]-50;g<-canvas(H)
    z<-res$contribution;z<-z[order(-z$Shapley_share),]
    header(g,'Two definitions of contribution',paste0('Aligned estimates with conditional 95% bootstrap intervals   /   ',
      format(BOOTSTRAP_B,big.mark=','),' resamples within Region'))
    g$line(95,76,110,76,ORANGE,.7);g$point(102.5,76,ORANGE,radius=2.5);g$text('Index variance share',116,76,9)
    g$line(240,76,255,76,BLUE,.7);g$point(247.5,76,BLUE,'square',2.3);g$text('Shapley R² allocation',261,76,9)
    left<-88;wid<-290;xmax<-43;x<-function(v)left+v/xmax*wid
    g$text('Estimate',410,88,8,MUTED,'centre');g$text('(%)',410,97,8,MUTED,'centre')
    for(i in 1:4) {
      yy<-129+(i-1)*42;g$line(left,yy,left+wid,yy,GRID,.5)
      g$text(z$dimension[i],left-4,yy,9.5,just='right')
      for(method in 1:2) {
        value<-100*if(method==1)z$covariance_share[i] else z$Shapley_share[i]
        low<-100*if(method==1)z$covariance_CI95_low[i] else z$Shapley_CI95_low[i]
        high<-100*if(method==1)z$covariance_CI95_high[i] else z$Shapley_CI95_high[i]
        dy<-if(method==1)-5.5 else 5.5;colour<-if(method==1)ORANGE else BLUE
        g$line(x(low),yy+dy,x(high),yy+dy,colour,.65)
        g$line(x(low),yy+dy-2.2,x(low),yy+dy+2.2,colour,.45)
        g$line(x(high),yy+dy-2.2,x(high),yy+dy+2.2,colour,.45)
        g$point(x(value),yy+dy,colour,if(method==1)'circle' else 'square',radius=1.95,border='white')
        g$text(sprintf('%.2f',value),401,yy+dy,8.6,just='centre')
      }
    }
    g$line(left,280,left+wid,280,RULE,.6)
    for(t in seq(0,40,10)){g$line(x(t),280,x(t),282.5,MUTED,.4);g$text(t,x(t),290,9,MUTED,'centre')}
    g$text('Contribution (%)',left+wid/2,311,9.5,just='centre')
    if(with_caption)g$caption(captions[4],351)
  }
  draw5<-function(with_caption=TRUE) {
    H<-if(with_caption)heights[5] else heights[5]-50;g<-canvas(H)
    header(g,'Prediction and regional transfer','Out-of-fold estimates from nested validation   /   separate leave-one-region-out sensitivity')
    g$text('a',65,64,10,bold=TRUE);g$text('Nested validation',79,64,10)
    tr<-axis(g,75,74,213,213,c(-1.65,1.65),c(-1.65,1.65),c(-1.5,-.5,.5,1.5),c(-1.5,-.5,.5,1.5),
      'Observed EMF average','Out-of-fold prediction',x_digits=1)
    g$line(tr$x(-1.65),tr$y(-1.65),tr$x(1.65),tr$y(1.65),'#A2AFBF',.7,lty='dashed')
    pp<-res$predictions;plot_points(g,pp$observed_EMF,pp$RF_OOF_prediction,pp$Region,tr)
    met<-res$RF_metrics
    g$text(sprintf('R² = %.3f',met$nested_OOF_R2),84,93,8.6)
    g$text(sprintf('RMSE = %.3f',met$nested_OOF_RMSE),84,104,8.6)
    g$text(sprintf('MAE = %.3f',met$nested_OOF_MAE),84,115,8.6)
    g$point(260,266,BLUE,radius=1.7,alpha=.65);g$text('GQ',270,266,8.5)
    g$point(260,278,ORANGE,'diamond',2.7,.8,'#B3822C');g$text('FT',270,278,8.5)
    g$text('b',355,83,10,bold=TRUE);g$text('Region transfer',369,83,10)
    tt<-res$transfer;tt<-tt[match(c('FT','GQ'),tt$test_region),]
    left<-356;wid<-121;xlim<-c(min(-.25,floor(min(tt$R2)*10)/10-.1),1)
    x<-function(v)left+(v-xlim[1])/diff(xlim)*wid
    g$line(x(0),104,x(0),212,RULE,.55,lty='dashed');g$line(left,212,left+wid,212,RULE,.6)
    for(t in c(0,.5,1)){g$line(x(t),212,x(t),215,MUTED,.4);g$text(sprintf('%.1f',t),x(t),222,9,MUTED,'centre')}
    for(i in 1:2) {
      yy<-131+(i-1)*64;g$text(paste0(tt$train_region[i],' to ',tt$test_region[i]),left-4,yy,8,just='right')
      g$point(x(tt$R2[i]),yy,if(i==1)BLUE else ORANGE,radius=2.25)
      g$text(sprintf('%.3f',tt$R2[i]),x(tt$R2[i]),yy-15,9,just='centre')
    }
    g$text('Transfer R²',left+wid/2,242,9.5,just='centre')
    g$text('Test-set RMSE',356,270,8.3,MUTED)
    for(i in 1:2)g$text(sprintf('%s to %s  %.3f',tt$train_region[i],tt$test_region[i],tt$RMSE[i]),356,282+(i-1)*12,8,MUTED)
    if(with_caption)g$caption(captions[5],372)
  }
  drawers<-list(draw1,draw2,draw3,draw4,draw5)
  stems<-c('Fig1_Correlation_heatmap','Fig2_Dimension_EMF_associations',
    'Fig3_RF_variable_importance','Fig4_Contribution_partitioning','Fig5_RF_validation')
  capdir<-file.path(out,'TIFF_带图注');puredir<-file.path(out,'TIFF_纯图');pngdir<-file.path(out,'PNG预览')
  for(dir in c(capdir,puredir,pngdir))dir.create(dir,showWarnings=FALSE)
  message('生成五幅原生R图：TIFF、PNG、可编辑矢量PPTX')
  for(i in 1:5) {
    grDevices::tiff(file.path(capdir,paste0(stems[i],'.tiff')),width=W/72,height=heights[i]/72,
      units='in',res=DPI,compression='lzw',type='cairo',bg='white')
    drawers[[i]](TRUE);grDevices::dev.off()
    grDevices::tiff(file.path(puredir,paste0(stems[i],'.tiff')),width=W/72,
      height=(heights[i]-if(i==3)52 else 50)/72,units='in',res=DPI,compression='lzw',type='cairo',bg='white')
    drawers[[i]](FALSE);grDevices::dev.off()
    grDevices::png(file.path(pngdir,paste0(stems[i],'.png')),width=W/72,height=heights[i]/72,
      units='in',res=200,type='cairo',bg='white')
    drawers[[i]](TRUE);grDevices::dev.off()
  }
  # 同参考PPTX使用13.333×10英寸画布，各图按自身比例放置，完整保留图注。
  template<-tempfile(fileext='.pptx');print(officer::read_pptx(),target=template)
  stage<-tempfile('EMF_ppt_template_');dir.create(stage);unzip(template,exdir=stage)
  xml<-xml2::read_xml(file.path(stage,'ppt','presentation.xml'))
  ns<-c(p='http://schemas.openxmlformats.org/presentationml/2006/main')
  sz<-xml2::xml_find_first(xml,'//p:sldSz',ns)
  xml2::xml_set_attr(sz,'cx','12192000');xml2::xml_set_attr(sz,'cy','9144000')
  xml2::write_xml(xml,file.path(stage,'ppt','presentation.xml'))
  zip::zipr(template,files=list.files(stage,recursive=TRUE,all.files=TRUE,no..=TRUE),root=stage,mode='mirror')
  ppt<-officer::read_pptx(template)
  for(i in 1:5) {
    ppt<-officer::add_slide(ppt,layout='Blank',master='Office Theme')
    width<-min(12.6,9.4*W/heights[i]);height<-width*heights[i]/W
    ppt<-officer::ph_with(ppt,rvg::dml(code=drawers[[i]](TRUE),bg='white',editable=TRUE),
      location=officer::ph_location(left=(13+1/3-width)/2,top=(10-height)/2,width=width,height=height))
    ppt<-officer::set_notes(ppt,value=c(captions[i],
      'R translation of the supplied Python statistical workflow. Numerical stochastic results are recomputed in R.',
      paste0('Data: ',INPUT_FILE),'https://imbs-hl.github.io/ranger/reference/ranger.html'),
      location=officer::notes_location_type('body'))
  }
  print(ppt,target=file.path(out,'EMF_visual_editable.pptx'))
  writeLines(unlist(lapply(captions,function(s)c(s,''))),file.path(out,'Figure_legends.txt'),useBytes=TRUE)
  writeLines(c('运行说明','',
    '在RStudio中以UTF-8打开EMF_statistics_visual.R并点击Source。脚本不调用Python。',
    paste0('默认输入：',INPUT_FILE),paste0('本次输出：',out),
    '与原Python一致：相关/BH校正、OLS和HC3、方差分解、Shapley、2000次区域内Bootstrap、',
    '一次5×4嵌套交叉验证、每折100次验证集置换、归一化MDI、跨区域验证和全局8NN残差诊断。',
    '第一次缺失包安装到任务专用EMF_R_library，后续使用已有包可离线运行。',
    '注意：R与Python的随机数、分折和随机森林后端不同，随机估计不保证逐位一致。',
    '图件使用同目录参考图的蓝橙配色、字体、点形、下三角和面板布局，数值由R重算。',
    '输出：5幅带图注TIFF、5幅纯图TIFF（600dpi LZW）、5幅PNG预览、5页可编辑矢量PPTX。',
    '完整统计对象为EMF_complete_results.rds；单文件汇总为EMF_statistics.txt；明细为统计明细CSV。',
    'EMF是其维度的构造型平均值，贡献与重要性不表示独立因果效应。',
    '颜色来源：Okabe–Ito https://jfly.uni-koeln.de/color/index.html',
    'RF API：https://imbs-hl.github.io/ranger/reference/ranger.html'),
    file.path(out,'README_中文.txt'),useBytes=TRUE)
  message('转换与出图完成：',normalizePath(out,winslash='/',mustWork=TRUE))
  invisible(list(output=out,results=res))
}
EMF_result <- main()
