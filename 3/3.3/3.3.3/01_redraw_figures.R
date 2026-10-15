# 仅更新图形，不重算模型、P值和自助区间。
# Rscript 01_redraw_figures.R "D:/Users/Lenovo/Desktop/原始/3.3/3.3.3"
options(stringsAsFactors=FALSE)
if(.Platform$OS.type=="windows") try(Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8"),silent=TRUE)
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else "D:/Users/Lenovo/Desktop/原始/3.3/3.3.3"
extra_lib <- Sys.getenv("VPA_R_LIB","")
if(nzchar(extra_lib)) .libPaths(c(extra_lib,.libPaths()))
if(dir.exists(file.path(out,"R_library"))) .libPaths(c(file.path(out,"R_library"),.libPaths()))
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
src <- tryCatch(sys.frame(1)$ofile,error=function(e) NULL)
ff <- grep("^--file=",commandArgs(),value=TRUE)
code_dir <- if(!is.null(src)) dirname(normalizePath(src)) else if(length(ff)) dirname(normalizePath(sub("^--file=","",ff[1]))) else getwd()
source(file.path(code_dir,"VPA_plot_functions.R"),local=TRUE,encoding="UTF-8")
candidates <- file.path(out,c("VPA_results.json","Supporting_results/VPA_results.json"))
found <- candidates[file.exists(candidates)]
if(!length(found)) stop("未找到分析JSON，请先运行 VPA_EMF_avg.R。")
r <- jsonlite::fromJSON(found[1],simplifyDataFrame=TRUE)
plots <- make_vpa_plots(r)
specs <- export_vpa_plots(plots,out)
if(requireNamespace("officer",quietly=TRUE)&&requireNamespace("rvg",quietly=TRUE)) {
  ppt <- officer::read_pptx()
  for(s in specs) {
    ppt <- officer::add_slide(ppt,layout="Blank",master="Office Theme")
    ppt <- officer::ph_with(ppt,rvg::dml(ggobj=s$plot),location=officer::ph_location(left=.1,top=.1,width=9.8,height=9.8*s$height/s$width))
  }
  print(ppt,target=file.path(out,"VPA_figures_R_editable.pptx"))
}
message("重绘完成：图1保持第一版；图2在第一版单面板和配色基础上改为横线长度对比。")
