
### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(4)
  library(ggplot2)

### load previously compiled data data
  load("~/dest2_glm_baypass_annotation.Rdata") ### made by `DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/baypass/collectBayPass.R`

### joint p-value threshold, no block bootstrap
  pvals <- expand.grid(sapply(c(-1:-6), function(x) c(1:9)*10^x))[,1]
  m[,xtx.rank:=rank(xtx.p)/(1+length(xtx.p))]
  m[,cont.rank:=rank(cont.p)/(1+length(cont.p))]

  enrich.stat <- foreach(p.i = pvals, .combine="rbind", .errorhandling="remove")%dopar%{
    foreach(p.j=pvals, .combine="rbind", .errorhandling="remove")%do%{
      #p.i<-.005; p.j<-.1
      message(paste(p.i, p.j, sep=" / "))
      xtx_cont.fet <- fisher.test(table(m$xtx.rank<=p.i, m$cont.rank<=p.j))
      data.table(xtx.thr=p.i, cont.thr=p.j, or=xtx_cont.fet$estimate, lci=xtx_cont.fet$conf.int[1], uci=xtx_cont.fet$conf.int[1])
    }
  }

  ggplot(data=enrich.stat[or!=0], aes(x=-log10(xtx.thr), y=-log10(cont.thr), fill=log2(or))) + geom_tile()

  m[xtx.p<1e-5 & cont.p<1e-5]
