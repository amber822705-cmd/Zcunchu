# EMF-nutrient 四组方差分解分析
# 运行: Rscript VPA_EMF_nutrient.R "D:/Users/Lenovo/Desktop/原始/new 模块.xlsx" "D:/Users/Lenovo/Desktop/原始/3.3/3.3.3"
# 可选第三参数 --analysis-only: 仅生成分析、CSV、JSON 和 TIFF。仅重绘请运行 01_redraw_figures.R。
# 配套：VPA_plot_functions.R，需与本脚本放在同一目录。
# R >= 4.5; 依赖采用 CRAN 官方包。所有随机过程均固定种子。
options(stringsAsFactors=FALSE, timeout=300)
if (.Platform$OS.type == "windows") try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"), silent=TRUE)
extra_lib <- Sys.getenv("VPA_R_LIB", "")
if (nzchar(extra_lib)) .libPaths(c(extra_lib, .libPaths()))
args <- commandArgs(trailingOnly=TRUE)
input <- if(length(args)>0) args[1] else "D:/Users/Lenovo/Desktop/原始/new 模块.xlsx"
out <- if(length(args)>1) args[2] else "D:/Users/Lenovo/Desktop/原始/3.3/VPA_EMF_nutrient"
analysis_only <- "--analysis-only" %in% args
local_lib <- file.path(out,"R_library")
if(dir.exists(local_lib)) .libPaths(c(local_lib,.libPaths()))
pkgs <- c("readxl","vegan","permute","ggplot2","patchwork","jsonlite")
if(!analysis_only) pkgs <- c(pkgs,"officer","rvg","flextable","openxlsx")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing_pkgs)) stop("请先运行 install.packages(c(",paste(sprintf('"%s"',missing_pkgs),collapse=","),"), repos='https://cloud.r-project.org')")
suppressPackageStartupMessages({library(vegan);library(ggplot2);library(patchwork)})
dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(20261014)
N_PERM <- 9999L
N_BOOT <- 2000L
cn <- c("Plot ID","Region","土壤含水率(%)","Salinity","pH","C:N比","N:P比","SOC","TN","TP","MAT_C（℃）","MAP_mm_yr（毫米 / 年）","Distance_coast_m","Distance_estuary_model_m","EMF-nutrient","EMF-nutrient","EMF-enzyme","EMF-vegetation","EMF-nutrient")
d0 <- as.data.frame(readxl::read_excel(input,sheet="计算",.name_repair="minimal"))
names(d0) <- trimws(names(d0))
if(anyDuplicated(names(d0))) stop("存在重复列名，请核查。")
miss <- setdiff(cn,names(d0)); if(length(miss)) stop("计算表缺列: ",paste(miss,collapse=", "))
d0 <- d0[,cn]
d0 <- d0[!is.na(d0[[1]]) & nzchar(trimws(as.character(d0[[1]]))),]
d0[[1]] <- trimws(as.character(d0[[1]])); d0[[2]] <- trimws(as.character(d0[[2]]))
if(anyDuplicated(d0[[1]])) stop("Plot ID 重复，请核实独立样本。")
for(j in 3:19) {
  old <- d0[[j]]; d0[[j]] <- suppressWarnings(as.numeric(old))
  if(any(!is.na(old) & is.na(d0[[j]]))) stop("非数值内容位于: ",cn[j])
}
need <- c(1:14,19)
keep <- complete.cases(d0[,need]) & apply(d0[,c(3:14,19)],1,function(z) all(is.finite(z)))
excluded <- d0[!keep,,drop=FALSE]; d <- d0[keep,,drop=FALSE]
if(nrow(d)<20) stop("完整样本过少，不能可靠拟合 12 个自变量。")
raw_x <- d[,3:14]; names(raw_x) <- c("Moisture","Salinity","pH","CN","NP","SOC","TN","TP","MAT","MAP","Coast","Estuary")
if(any(vapply(raw_x,sd,numeric(1))==0)) stop("存在常数自变量，请先核查。")
y <- d[["EMF-nutrient"]]; if(sd(y)==0) stop("响应变量无变异。")
region <- factor(d$Region)
groups <- list(A=1:3,B=4:5,C=6:8,D=9:12)
group_names <- c(A="土壤生境组",B="化学计量组",C="植物组",D="区域环境组")
colors <- c(A="#0072B2",B="#D55E00",C="#009E73",D="#CC79A7")
X <- scale(raw_x)
bits <- function(m) which(as.logical(intToBits(m)[1:4]))
ord <- (1:15)[order(vapply(1:15,function(m) length(bits(m)),integer(1)),1:15)]
labels <- vapply(ord,function(m) paste(LETTERS[bits(m)],collapse="∩"),character(1))
adj_lm <- function(yy,xx,z=NULL) {
  mm <- cbind(Intercept=1,z,xx); fit <- lm.fit(mm,yy); n <- length(yy)
  r2 <- 1-sum(fit$residuals^2)/sum((yy-mean(yy))^2)
  if(n<=fit$rank) stop("残差自由度不足。")
  c(R2=r2,Adj_R2=1-(1-r2)*(n-1)/(n-fit$rank),Df=fit$rank-1)
}
partition <- function(yy,xx,z=NULL) {
  fits <- t(vapply(1:15,function(m) adj_lm(yy,xx[,unlist(groups[bits(m)]),drop=FALSE],z),numeric(3)))
  base <- if(is.null(z)) 0 else adj_lm(yy,matrix(numeric(0),nrow=length(yy),ncol=0),z)["Adj_R2"]
  rr <- c(0,fits[,"Adj_R2"]-base)
  frac <- numeric(15)
  for(m in ord) {
    proper <- (1:15)[vapply(1:15,function(k) bitwAnd(k,m)==k && k!=m,logical(1))]
    frac[m] <- rr[16]-rr[bitwXor(15L,m)+1]-sum(frac[proper])
  }
  stopifnot(abs(sum(frac)-rr[16])<1e-10)
  list(fits=fits,frac=frac,total=rr[16],base=base)
}
pa <- partition(y,X)
vp <- vegan::varpart(y,as.data.frame(X[,groups$A]),as.data.frame(X[,groups$B]),as.data.frame(X[,groups$C]),as.data.frame(X[,groups$D]))
# 独立检验: 根据 vegan 的前四个基础模型自动识别 a-o 的组成员，逐项核对全部 15 分量。
vletters <- letters[1:15]
memletters <- lapply(rownames(vp$part$fract)[1:4],function(s) strsplit(sub(".*\\[([a-o]+)\\].*","\\1",s),"")[[1]])
vmasks <- vapply(vletters,function(l) sum(2^(which(vapply(memletters,function(z) l %in% z,logical(1)))-1)),numeric(1))
vcol <- grep("^Adj\\.R\\.squar",names(vp$part$indfract),value=TRUE)
if(length(vcol)!=1L) stop("无法识别vegan校正R²列。")
vf <- vp$part$indfract[1:15,vcol]
stopifnot(length(vf)==15L,all(is.finite(vf)),max(abs(pa$frac[vmasks]-vf))<1e-8)
if(max(abs(pa$fits[,"Adj_R2"]-vapply(1:15,function(m) vegan::RsquareAdj(vegan::rda(y,X[,unlist(groups[bits(m)]),drop=FALSE]))$adj.r.squared,numeric(1))))>1e-8) stop("R² 复核失败。")
perms_block <- permute::shuffleSet(nrow(d),nset=N_PERM,control=permute::how(blocks=region,nperm=N_PERM))
perms_free <- permute::shuffleSet(nrow(d),nset=N_PERM,control=permute::how(nperm=N_PERM))
test_rda <- function(xx,zz=NULL,pm=perms_block) {
  fit <- if(is.null(zz)) vegan::rda(y,xx) else vegan::rda(y,xx,zz)
  a <- anova(fit,permutations=pm,model="reduced")
  c(Df=a$Df[1],F=a$F[1],P=a$`Pr(>F)`[1])
}
fulltest <- test_rda(X)
puretests <- t(vapply(seq_along(groups),function(g) test_rda(X[,groups[[g]],drop=FALSE],X[,-groups[[g]],drop=FALSE]),numeric(3)))
freetests <- t(vapply(seq_along(groups),function(g) test_rda(X[,groups[[g]],drop=FALSE],X[,-groups[[g]],drop=FALSE],perms_free),numeric(3)))
pure <- data.frame(Group=names(groups),Group_name=unname(group_names),N_variables=lengths(groups),Adj_R2=pa$frac[c(1,2,4,8)],Percent=100*pa$frac[c(1,2,4,8)],Df=puretests[,"Df"],F=puretests[,"F"],P_block=puretests[,"P"],Q_BH=p.adjust(puretests[,"P"],method="BH"),P_free=freetests[,"P"])
# 在 Region 内进行分层自助抽样。该区间以样点条件独立为前提，不校正未知的空间自相关。
set.seed(20261015)
strata <- split(seq_len(nrow(d)),region)
boot <- replicate(N_BOOT,{
  ii <- unlist(lapply(strata,function(j) sample(j,length(j),replace=TRUE)),use.names=FALSE)
  b <- partition(y[ii],X[ii,,drop=FALSE]); c(b$frac[c(1,2,4,8)],b$total,sum(b$frac[-c(1,2,4,8)]))
})
ci <- t(apply(boot,1,quantile,probs=c(.025,.975),na.rm=TRUE))
pure$CI_low <- ci[1:4,1]; pure$CI_high <- ci[1:4,2]
fractions <- data.frame(Code=vletters[match(ord,vmasks)],Membership=labels,Order=vapply(ord,function(m) length(bits(m)),integer(1)),Adj_R2=pa$frac[ord],Percent=100*pa$frac[ord],Testable=vapply(ord,function(m) length(bits(m))==1,logical(1)))
for(i in 1:4) fractions[[LETTERS[i]]] <- vapply(ord,function(m) i %in% bits(m),logical(1))
unions <- data.frame(Membership=vapply(1:15,function(m) paste(LETTERS[bits(m)],collapse="+"),character(1)),pa$fits,check.names=FALSE)
model_summary <- data.frame(Item=c("样本数","自变量数","模型自由度","残差自由度","原始R2","校正R2","独立解释度合计","共同解释度合计","未解释比例","整体F","区域内置换P"),Value=c(nrow(d),ncol(X),pa$fits[15,"Df"],nrow(d)-pa$fits[15,"Df"]-1,pa$fits[15,"R2"],pa$total,sum(pure$Adj_R2),sum(fractions$Adj_R2[fractions$Order>1]),1-pa$total,fulltest["F"],fulltest["P"]))
# 敏感性1: 对 N:P 取自然对数，对两种距离取 log1p，再重新标准化。保留其余变量及样本。
if(any(raw_x$NP<=0) || any(as.matrix(raw_x[,c("Coast","Estuary")])<0)) stop("敏感性转换定义域不满足。")
lx <- raw_x; lx$NP <- log(lx$NP); lx$Coast <- log1p(lx$Coast); lx$Estuary <- log1p(lx$Estuary); LX <- scale(lx)
plog <- partition(y,LX)
logtests <- t(vapply(seq_along(groups),function(g) test_rda(LX[,groups[[g]],drop=FALSE],LX[,-groups[[g]],drop=FALSE]),numeric(3)))
# 敏感性2: 将 Region 作为条件协变量。独立分量定义为含 Region 的全模型与删组模型校正R²之差，分母仍为原始响应总方差。
Zregion <- model.matrix(~region)[,-1,drop=FALSE]
preg <- partition(y,X,Zregion)
regtests <- t(vapply(seq_along(groups),function(g) test_rda(X[,groups[[g]],drop=FALSE],cbind(Zregion,X[,-groups[[g]],drop=FALSE])),numeric(3)))
sensitivity <- rbind(
  data.frame(Analysis="主分析",Group=names(groups),Adj_R2=pure$Adj_R2,P_block=pure$P_block,Q_BH=pure$Q_BH),
  data.frame(Analysis="对数转换",Group=names(groups),Adj_R2=plog$frac[c(1,2,4,8)],P_block=logtests[,"P"],Q_BH=p.adjust(logtests[,"P"],"BH")),
  data.frame(Analysis="控制Region",Group=names(groups),Adj_R2=preg$frac[c(1,2,4,8)],P_block=regtests[,"P"],Q_BH=p.adjust(regtests[,"P"],"BH")))
sensitivity$Percent <- 100*sensitivity$Adj_R2
vif <- as.numeric(vegan::vif.cca(vegan::rda(y,X)))
descriptive <- do.call(rbind,lapply(c(3:14,19),function(j) {z<-d[[j]];data.frame(Variable=cn[j],N=length(z),Missing=sum(is.na(d0[[j]])),Mean=mean(z),SD=sd(z),Min=min(z),Median=median(z),Max=max(z),Skew=mean(((z-mean(z))/sqrt(mean((z-mean(z))^2)))^3))}))
diagnostics <- data.frame(Variable=cn[3:14],Short_name=colnames(X),Group=rep(names(groups),lengths(groups)),VIF=vif)
lmfull <- lm(y~X); cooks <- cooks.distance(lmfull); h <- hatvalues(lmfull)
obs <- data.frame(Plot_ID=d[[1]],Region=d$Region,Observed=y,Fitted=fitted(lmfull),Residual=resid(lmfull),Standardized_residual=rstandard(lmfull),Cook_D=cooks,Leverage=h,Influence_flag=cooks>4/nrow(d))
# 影响点敏感性: 仅用于评估结论稳定性，主分析不删除观测。
ii <- which(cooks<=4/nrow(d)); pinfl <- partition(y[ii],X[ii,,drop=FALSE])
influence_sensitivity <- data.frame(Group=names(groups),Main=pure$Adj_R2,Without_influential=pinfl$frac[c(1,2,4,8)],N_without=length(ii))
bp <- lm(resid(lmfull)^2~X); bp_stat <- nrow(d)*summary(bp)$r.squared
assumptions <- data.frame(Check=c("Shapiro-Wilk residual P","Breusch-Pagan auxiliary P","Max Cook D","Cook D > 4/n count","Max VIF","EMF mean check max error"),Value=c(shapiro.test(resid(lmfull))$p.value,pchisq(bp_stat,ncol(X),lower.tail=FALSE),max(cooks),sum(cooks>4/nrow(d)),max(vif),max(abs(d[[19]]-rowMeans(d[,15:18])),na.rm=TRUE)))
tables <- list(Model_summary=model_summary,Pure_effects=pure,All_fractions=fractions,Union_models=unions,Sensitivity=sensitivity,Influence_sensitivity=influence_sensitivity,Descriptive=descriptive,Collinearity=diagnostics,Assumptions=assumptions,Sample_diagnostics=obs,Analysis_data=d,Excluded_data=excluded)
write_utf8_csv <- function(z,p) {con<-file(p,open="wb");writeBin(charToRaw("\ufeff"),con);close(con);con<-file(p,open="a",encoding="UTF-8");write.table(z,con,sep=",",row.names=FALSE,col.names=TRUE,na="",qmethod="double");close(con)}
dir.create(file.path(out,"CSV"),showWarnings=FALSE)
for(nm in names(tables)) write_utf8_csv(tables[[nm]],file.path(out,"CSV",paste0(nm,".csv")))
correlation <- cor(X)
write_utf8_csv(data.frame(Variable=cn[3:14],correlation,check.names=FALSE),file.path(out,"CSV","Correlation.csv"))
res <- list(source=normalizePath(input,winslash="/"),sheet="计算",source_md5=unname(tools::md5sum(input)),n=nrow(d),excluded=nrow(excluded),region_counts=as.list(table(region)),nperm=N_PERM,nboot=N_BOOT,groups=lapply(seq_along(groups),function(i) list(code=names(groups)[i],name=group_names[i],variables=cn[groups[[i]]+2])),model=tables$Model_summary,pure=pure,fractions=fractions,unions=unions,sensitivity=sensitivity,influence_sensitivity=influence_sensitivity,descriptive=descriptive,collinearity=diagnostics,assumptions=assumptions,observations=obs,correlation=unname(correlation),total_adj=pa$total,total_r2=pa$fits[15,"R2"],unique_sum=sum(pure$Adj_R2),shared_sum=sum(fractions$Adj_R2[fractions$Order>1]),residual=1-pa$total,full_F=unname(fulltest["F"]),full_P=unname(fulltest["P"]),log_total=plog$total,region_base=preg$base,region_added=preg$total,verified_against_vegan=TRUE)
jsonlite::write_json(res,file.path(out,"VPA_results.json"),pretty=TRUE,auto_unbox=TRUE,digits=15,na="null")
saveRDS(list(varpart=vp,lm=lmfull,tables=tables,results=res,bootstrap=boot),file.path(out,"VPA_models.rds"))
capture.output(sessionInfo(),file=file.path(out,"R_sessionInfo.txt"))
capture.output(print(vp),file=file.path(out,"VPA_vegan_output.txt"))

# 图1沿用原版，图2使用优化版；两图共用绘图函数，保持统计值不变。
src <- tryCatch(sys.frame(1)$ofile,error=function(e) NULL)
ff <- grep("^--file=",commandArgs(),value=TRUE)
code_dir <- if(!is.null(src)) dirname(normalizePath(src)) else if(length(ff)) dirname(normalizePath(sub("^--file=","",ff[1]))) else getwd()
plot_helper <- file.path(code_dir,"VPA_plot_functions.R")
if(!file.exists(plot_helper)) stop("请将 VPA_plot_functions.R 与本脚本放在同一目录。")
source(plot_helper,local=TRUE,encoding="UTF-8")
plots <- make_vpa_plots(res)
plot_specs <- export_vpa_plots(plots,out)
mainplot <- plots$main
pS <- plots$sensitivity

notes <- c(
 "图1恢复原版：四组独立解释度及完整95%自助区间，配合15个有符号精确分量和成员点阵。图1的轴范围、配色及标注沿用原版。",
 "图2在第一版单面板基础上，将点改为从0延伸至估计值的横线；四组纵轴及按分析方案区分的配色保持第一版设置，横线长度用于比较各方案解释度。",
 "方法：EMF-nutrient 为单变量连续响应，使用多元线性回归/偏回归意义下的 VPA；不对响应进行 Hellinger 或对数转换。主分析12个自变量仅作z标准化，保留原始线性关系。",
 "A=土壤生境组（含水率、Salinity、pH），B=化学计量组（C:N、N:P），C=植物组（SOC、TN、TP），D=区域环境组（MAT、MAP、两种距离）。组名及成员按指定方案保留。",
 "SOC、TN、TP的采样组织和单位需结合实验记录核实。若其为土壤指标，植物组应更名为土壤碳氮磷组；若其与EMF构建使用相同原始指标，还需评估数学耦合，不能直接解读为因果驱动。",
 "独立解释度是全模型与删去该组模型的校正R²差。共同分量通过全部15个组联合模型的容斥关系计算，逐项与vegan::varpart交叉核验。",
 "独立与共同分量均相对于EMF-nutrient总变异。总解释率=15分量之和；未解释率=1-全模型校正R²。各组单独模型解释率有重叠，不能相加。",
 "负分量原样保留，不截零、不重新归一化。共同分量是共享统计解释，不能直接作为生态交互作用或协同作用证据。",
 "主检验为Region内限制置换（9999次），偏模型置换条件模型残差（model='reduced'）。四个独立分量的P值作BH校正；整体模型P单独报告。共享分量不能直接单独检验，不提供伪造的P值。",
 "95%区间由Region内按样点有放回分层抽样2000次得到，为条件于现有区域的描述性不确定性，不覆盖区域外推广误差；仅有两个Region，不能可靠估计跨区域随机效应。",
 "数据无坐标和更细的独立采样单元标识，区域限制置换不能消除全部空间自相关。P值和自助区间仍以区域内样点的交换性/条件独立为前提。",
 "土壤含水率存在超过100%的值。干质量基准含水率可以超过100%；因未提供计量基准，本分析保留原值，不擅自裁剪或改单位。",
 "敏感性分析的对数转换、Region条件化及Cook距离筛查均不替代主分析。Cook D>4/n只是影响点筛查阈值，不是删除样本的依据。",
 "参考：Peres-Neto et al. (2006), Ecology 87:2614–2625, doi:10.1890/0012-9658(2006)87[2614:VPSDMA]2.0.CO;2。",
 "软件说明：https://vegandevs.github.io/vegan/reference/varpart.html；https://vegandevs.github.io/vegan/reference/anova.cca.html；https://vegandevs.github.io/vegan/reference/RsquareAdj.html"
)
writeLines(enc2utf8(notes),file.path(out,"方法与图表说明.txt"),useBytes=TRUE)
if(!analysis_only) {
  # 独立重跑时，直接由R生成Excel、可编辑矢量PPTX和Word三线表。
  wb <- openxlsx::createWorkbook()
  for(nm in names(tables)) {
    openxlsx::addWorksheet(wb,nm,gridLines=FALSE)
    openxlsx::writeData(wb,nm,tables[[nm]],keepNA=FALSE)
    openxlsx::freezePane(wb,nm,firstRow=TRUE)
    openxlsx::setColWidths(wb,nm,cols=seq_len(ncol(tables[[nm]])),widths="auto")
    header <- openxlsx::createStyle(fontName="Microsoft YaHei",textDecoration="bold",border=c("top","bottom"),borderStyle=c("medium","thin"))
    bottom <- openxlsx::createStyle(border="bottom",borderStyle="medium")
    openxlsx::addStyle(wb,nm,header,rows=1,cols=seq_len(ncol(tables[[nm]])),gridExpand=TRUE)
    if(nrow(tables[[nm]])) openxlsx::addStyle(wb,nm,bottom,rows=nrow(tables[[nm]])+1,cols=seq_len(ncol(tables[[nm]])),gridExpand=TRUE)
  }
  openxlsx::addWorksheet(wb,"Notes");openxlsx::writeData(wb,"Notes",data.frame(Note=notes))
  openxlsx::saveWorkbook(wb,file.path(out,"VPA_results_R.xlsx"),overwrite=TRUE)
  ppt <- officer::read_pptx()
  for(s in plot_specs) {
    ppt <- officer::add_slide(ppt,layout="Blank",master="Office Theme")
    ppt <- officer::ph_with(ppt,rvg::dml(ggobj=s$plot),location=officer::ph_location(left=.1,top=.1,width=9.8,height=9.8*s$height/s$width))
  }
  print(ppt,target=file.path(out,"VPA_figures_R_editable.pptx"))
  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc,"EMF-nutrient 方差分解分析",style="heading 1")
  doc <- officer::body_add_par(doc,sprintf("完整样本n=%d，全模型校正R²=%.2f%%，独立分量合计%.2f%%，共同分量合计%.2f%%。",nrow(d),100*pa$total,100*res$unique_sum,100*res$shared_sum))
  report_tabs <- list("表1 四组独立解释度"=pure[,c("Group","Group_name","Percent","F","P_block","Q_BH","CI_low","CI_high")],"表2 全部独立与共同分量"=fractions[,1:6],"表3 敏感性分析"=sensitivity,"表4 自变量描述统计"=descriptive,"表5 共线性诊断"=diagnostics)
  for(nm in names(report_tabs)) {
    doc<-officer::body_add_par(doc,nm,style="heading 2")
    ft<-flextable::flextable(report_tabs[[nm]])
    ft<-flextable::theme_booktabs(ft);ft<-flextable::font(ft,fontname="Microsoft YaHei",part="all");ft<-flextable::fontsize(ft,size=9,part="all");ft<-flextable::autofit(ft)
    ft<-flextable::set_table_properties(ft,layout="autofit",width=1)
    doc<-flextable::body_add_flextable(doc,ft)
    doc<-officer::body_add_par(doc,"注：解释度以总变异为分母。独立分量P值为区域内置换，Q为四组BH校正；共同分量不单独检验。CI为R²原比例。")
  }
  for(s in c(notes,plots$captions)) doc<-officer::body_add_par(doc,s)
  print(doc,target=file.path(out,"VPA_tables_R.docx"))
}
cat(sprintf("VPA完成: n=%d, R2=%.6f, adjusted R2=%.6f\n",nrow(d),pa$fits[15,"R2"],pa$total))
print(pure,digits=5)
cat("共同分量合计 =",res$shared_sum,"；未解释 =",res$residual,"\n")
