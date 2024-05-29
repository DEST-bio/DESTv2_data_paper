library(data.table)
library(ggplot2)
library(patchwork)
load("/Users/alanbergland/dest2_seasonality/C2_windowPerms_pod.Rdata")

nSNP_thr=150

C2.wza.plot <- ggplot(data=win.out[nSNPs>nSNP_thr][perm==0]) +
geom_line(data=win.out[perm==0], aes(x=pos_mean/1e6, y=-C2.wZa.p), color="black") +
facet_grid(~chr, scales="free_x") + ylab("-log10(C2 wZa p)") + xlab("Pos (Mb)") + theme_bw()

C2.wza.pod.plot <- ggplot(data=win.out[nSNPs>nSNP_thr][perm==0]) +
geom_line(data=win.out[perm==0],aes(x=pos_mean/1e6, y=-C2.wZa.pod.p), color="black") +
facet_grid(~chr, scales="free_x") + ylab("-log10(C2 wZa p pod)") + xlab("Pos (Mb)") + theme_bw() +
geom_hline(yintercept=-log10(.0001/dim(win.out)[1]))

C2.wza.pod.orig.plot <- ggplot(data=win.out[nSNPs>nSNP_thr][perm==0]) +
geom_line(data=win.out[perm==0],aes(x=pos_mean/1e6, y=-C2.wZa.pod.orig.p), color="black") +
facet_grid(~chr, scales="free_x") + ylab("-log10(C2 wZa p pod orig)") + xlab("Pos (Mb)") + theme_bw() +
geom_hline(yintercept=-log10(.0001/dim(win.out)[1]))

C2.wza.plot /C2.wza.pod.plot /C2.wza.pod.orig.plot
ggplot(data=win.out, aes(x=C2.wZa.pod, y=C2.wZa)) + geom_point()
