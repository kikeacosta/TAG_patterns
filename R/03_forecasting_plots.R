if(mg<=25){
  p <- ggplot(DF, aes(x=years1, y=eta1))+
    geom_line(data=filter(DF, type=="Fitted+Forecast"),
              aes(y=eta1, col=colo, group=ages),size=1)+
    geom_point(data=filter(DF, type=="Actual grouped"),
               aes(y=eta1, col=colo))+
    facet_wrap(~group, 5, 5, scales="free_y")+
    theme(legend.position = "none")+
    scale_color_viridis(discrete=FALSE, option="inferno")+
    labs(x="year", y="log-mortality", title=cou.j)+
    geom_rect(data=DF[1:mg,],
              aes(xmin = min(t2)-0.2, xmax = max(t2)+0.2,
                  ymin = -Inf, ymax = Inf),
              fill="cyan", alpha=0.05)
  p
}else{
  DFsel <- subset(DF, ages%in%seq(0,90,10))
  p <- ggplot(DFsel, aes(x=years1, y=eta1))+
    geom_line(data=filter(DFsel, type=="Fitted+Forecast"),
              aes(y=eta1, col=colo, group=ages),size=1)+
    geom_point(data=filter(DFsel,type=="Actual grouped"),
               aes(y=eta1, col=colo))+
    facet_wrap(~group, 2, 5, scales="free_y")+
    theme(legend.position = "none")+
    scale_color_viridis(discrete=FALSE, option="inferno")+
    labs(x="year", y="log-mortality", title=cou.j)+
    geom_rect(data=DFsel[1:mg,],
              aes(xmin = min(t2)-0.2, xmax = max(t2)+0.2,
                  ymin = -Inf, ymax = Inf),
              fill="cyan", alpha=0.05)
  p
}


p <- ggplot(DF, aes(x=ages, y=eta1, color=type)) +
  geom_segment(data=filter(DF, type=="Actual grouped"),
               aes(x=ages, y=eta1, xend=ages.up, yend=eta1.up), linewidth=1)+
  geom_line(data=filter(DF, type=="Fitted+Forecast"),
            aes(y=eta1),linewidth=1)+
  geom_line(data=filter(DF, type=="Fitted"),
            aes(y=eta1),linewidth=1, linetype="dotted")+
  facet_wrap(~years1, 2, 6, scales="free_y")+
  labs(x="age", y="log-mortality", title=cou.j)


# check to see if repeated
if(mg<=25){
  p <- ggplot(DF, aes(x=years1, y=eta1)) +
    geom_line(data=filter(DF, type=="Fitted"),
              aes(y=eta1, col=colo, group=ages),size=1)+
    geom_point(data=filter(DF, type=="Actual grouped"),
               aes(y=eta1, col=colo))+
    facet_wrap(~group, 5, 5, scales="free_y")+
    theme(legend.position = "none") +
    scale_color_viridis(discrete=FALSE, option="inferno")+
    labs(x="year", y="log-mortality", title=cou.j)
  p
}else{
  DFsel <- subset(DF, ages%in%seq(0,90,10))
  p <- ggplot(DFsel, aes(x=years1, y=eta1)) +
    geom_line(data=filter(DFsel, type=="Fitted"),
              aes(y=eta1, col=colo, group=ages),size=1)+
    geom_point(data=filter(DFsel, type=="Actual grouped"),
               aes(y=eta1, col=colo))+
    facet_wrap(~group, 2, 5, scales="free_y")+
    theme(legend.position = "none") +
    scale_color_viridis(discrete=FALSE, option="inferno")+
    labs(x="year", y="log-mortality", title=cou.j)
  p
}