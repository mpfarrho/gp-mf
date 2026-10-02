# ---------------------------------------------------------------------------------------------
# example: nowcast and one-quarter-ahead forecast of US real GDP growth (GDPC1) and GDP deflator
# inflation (GDPCTPI) with the small information set, using data pulled from ALFRED
#
# run with the folder of this script as working directory
#
# differences to the data used in the paper (FRED-MD/QD):
#  - the S&P 500 is not included (ALFRED only provides the last ten years), so K = 11
#  - initial claims (ICSA, CLAIMSx in FRED-MD) start in 1967, which determines the sample start
#
# timing: the nowcast refers to the quarter of the vintage date, the forecast to the quarter after
# that; a backcast of the previous quarter is added in case its GDP figures have not been released
# yet; all use the data released until the vintage date
library(alfred)
library(zoo)
library(lubridate)
library(ggplot2)

source("aux_func.R")
source("gpmf_func.R")
source("facpc_func.R")

set.seed(1)

# ---------------------------------------------------------------------------------------------
# settings
vintage <- Sys.Date()  # data as available on this day (for a past vintage use, e.g., as.Date("2024-05-15"))

run.type <- "xalm"     # mixed-frequency compression: "u", "f", "avg", "br", "xalm", "alm_3", "leg_3", "leg_5", "ber_3", "ber_5", "fou_3"
run.model <- "gp-sqexp" # conditional mean: "gp-sqexp", "bart-250", "lin-hs"
run.sv <- TRUE         # stochastic volatility

sl.target <- c("GDPC1","GDPCTPI")

# set lag structure of MIDAS design matrix
p.y <- 4
p.x <- 4
cons <- FALSE # demeaned variables

# short run for illustration (paper: nburn <- 3000; nsave <- 3000; nthin <- 3)
nburn <- 1000
nsave <- 1000
nthin <- 1

tau.smry <- c(0.05, 0.1, 0.16, 0.5, 0.84, 0.9, 0.95)

run.mod <- unlist(strsplit(run.model, "-")); run.spec <- run.mod[[2]]; run.mod <- run.mod[[1]]

midas.par1 <- midas.par2 <- 0
par1_pr <- par2_pr <- 0

shrink <- TRUE
if(run.mod == "bart"){
  n.trees <- as.numeric(run.spec)
  gp.type <- ""
}else if(run.mod == "gp"){
  n.trees <- NA
  gp.type <- run.spec
}else if(run.mod == "lin"){
  n.trees <- NA
  gp.type <- ""
  if(run.spec == "flat"){
    shrink <- FALSE
  }else if(run.spec == "hs"){
    shrink <- TRUE
  }
}

# ---------------------------------------------------------------------------------------------
# download data and transform as in FRED-MD
# (1) no transformation; (2) diff; (4) log; (5) 100 x diff of log; (6) diff of 100 x diff of log
sl.vars <- c("INDPRO" = 5, "DPCERA3M086SBEA" = 5,                             # production/real
             "PAYEMS" = 5, "UNRATE" = 2, "ICSA" = 5, "CES0600000007" = 1,     # labour market
             "HOUST" = 4,                                                     # housing
             "CPIAUCSL" = 5,                                                  # prices
             "FEDFUNDS" = 2, "GS10" = 2, "BAAFFM" = 1                         # financial markets
)

get_vintage <- function(id, freq){
  # missing values in ALFRED trigger a coercion warning
  x <- suppressWarnings(get_alfred_series(id, "value", observation_start = "1959-01-01",
                                          realtime_start = as.character(vintage), realtime_end = as.character(vintage)))
  if(id == "ICSA"){ # weekly series: averages of completed months
    x <- x[x$date < floor_date(vintage, "month"),]
    x <- aggregate(value ~ date, data = data.frame("date" = floor_date(x$date, "month"), "value" = x$value), FUN = mean)
  }
  x$value <- approx(x$date, x$value, xout = x$date)$y # linear interpolation of missings within the sample
  ts(x$value, start = c(year(x$date[1]), (month(x$date[1]) - 1) / (12/freq) + 1), frequency = freq)
}

data_raw <- list()
for(j in names(sl.vars)){
  y <- get_vintage(j, freq = 12)
  if(sl.vars[[j]]==2){
    y <- diff(y)
  }else if(sl.vars[[j]]==4){
    y <- log(y)
  }else if(sl.vars[[j]]==5){
    y <- 100*diff(log(y))
  }else if(sl.vars[[j]]==6){
    y <- diff(100*diff(log(y)))
  }
  data_raw[[j]] <- y
}
for(j in sl.target){
  y <- get_vintage(j, freq = 4)
  data_raw[[j]] <- 100*(((y/stats::lag(y,-1))^4)-1) # compounded annual rate of change
}

# balanced panels: monthly predictors until the last month with all series released, targets until
# the last released quarter
Xraw <- do.call(ts.intersect,data_raw[names(sl.vars)])
Yraw <- do.call(ts.intersect,data_raw[sl.target])
colnames(Xraw) <- names(sl.vars)
N <- ncol(Yraw)
K <- ncol(Xraw)
Tq <- nrow(Yraw)
Tm <- nrow(Xraw)

# ---------------------------------------------------------------------------------------------
# real-time calendar: the nowcast refers to the quarter of the vintage date, the forecast to the
# quarter after that and the backcast to the quarter before
run.ho <- as.numeric(as.yearqtr(vintage))
message("Vintage: ", vintage, ". Quarterly data until ", as.yearqtr(tsp(Yraw)[2] + 1e-8), ", monthly data until ", as.yearmon(tsp(Xraw)[2] + 1e-8), ".")

# hh = 0 nowcast, hh = 1 forecast, hh = -1 backcast (in case the previous quarter has not been released yet)
grid.hh <- c(if(round(4*(run.ho - tsp(Yraw)[2])) > 1) -1, 0, 1)

fcst.out <- fcst.draws <- list()
for(hh in grid.hh){
  # ---------------------------------------------------------------------------------
  # direct predictive regressions: the publication lags of the target quarter are imposed on
  # each quarter of the sample
  lag.q <- round(4*(run.ho + hh/4 - tsp(Yraw)[2]))         # quarters between target quarter and last release of the targets
  lag.m <- round(12*(run.ho + hh/4 + 2/12 - tsp(Xraw)[2])) # months between final month of target quarter and last monthly observation

  # row of the most recent monthly observation available for each quarter
  id.m <- round(12*(as.numeric(time(Yraw)) + 2/12 - tsp(Xraw)[1])) + 1 - lag.m
  sl.t <- which((1:Tq) - lag.q - (p.y-1) >= 1 & id.m - (3*p.x-1) >= 1)

  # lags of the targets (by lag) and the 3*p.x most recent monthly observations (by predictor)
  Ylags <- t(sapply(sl.t, function(tt) c(t(Yraw[tt - lag.q - (0:(p.y-1)),]))))
  Yth <- c(t(Yraw[Tq - (0:(p.y-1)),]))
  colnames(Ylags) <- names(Yth) <- paste0(rep(sl.target, p.y),"_",rep(1:p.y, each = N))

  X_m <- t(sapply(sl.t, function(tt) c(Xraw[id.m[tt] - (0:(3*p.x-1)),])))
  Xth_m <- c(Xraw[Tm - (0:(3*p.x-1)),])
  colnames(X_m) <- names(Xth_m) <- paste0(rep(colnames(Xraw), each = 3*p.x),"_",1:(3*p.x))

  # compression of MIDAS design matrix with factor structure
  if(run.type == "f"){
    X_tmp <- rbind(X_m, Xth_m)

    fnum <- try(as.numeric(get_num_factors(scale(X_tmp), kmax = 20, criteria = "IC2")$ic))
    if(is(fnum,'try-error')) fnum <- 4
    fnum[fnum > 20] <- 20; fnum[fnum < 3] <- 3

    # draw factors with PCA
    X_facs <- try(get_em_factors(X_tmp, n = fnum, it_max = 50)$factors, silent = TRUE)
    if(is(X_facs,'try-error')){
      X_facs <- matrix(scale(prcomp(X_tmp, center = TRUE, scale = TRUE)$x[,1:fnum]), nrow(X_tmp), ncol = fnum)
    }
    colnames(X_facs) <- paste0("FAC",1:fnum)

    X_m <- X_facs[1:nrow(X_m),]
    Xth_m <- X_facs[nrow(X_facs),]
  }

  # final design matrices
  y <- Yraw[sl.t,,drop=FALSE]
  X <- cbind(Ylags,X_m)
  Xth <- matrix(c(Yth,Xth_m),nrow=1); colnames(Xth) <- colnames(X)

  # ---------------------------------------------------------------------------------------------
  # estimation
  prior_setup <- list(set.mean=run.mod,set.sv=run.sv,set.midas=run.type,
                      p.y=p.y,K=K,p.x=p.x,cons=cons,

                      a0_hom=3,b0_hom=0.3,sv_sig=0.1, # variances
                      shrink=shrink,midas.par1=midas.par1,midas.par2=midas.par2,par1_pr=par1_pr,par2_pr=par2_pr, # linear priors
                      cgm.level=0.95,cgm.exp=2,sd.mu=1.96,num.trees=n.trees, # bart priors

                      gp.type = gp.type, xi=1,h= 0.1,zeta=7/2, zeta.grid = seq(3/2, 15/2, 1),
                      max.try.solve=10, samp.gp.hyper=TRUE # gp setup
  )
  sim <- npmf(y=y,X=X,Xth=Xth,prior_setup=prior_setup,nburn=nburn,nsave=nsave,nsave_lim=min(nsave,1000),nthin=nthin)

  # predictive distribution (annualized growth rates in percent)
  fcst <- sim$fcst
  fcst.smry <- rbind(apply(fcst, 2, quantile, tau.smry, na.rm = TRUE),
                     "sd"= apply(fcst, 2, sd, na.rm = TRUE))
  colnames(fcst.smry) <- sl.target

  fcst.label <- paste0(c("backcast ","nowcast ","forecast ")[hh+2], as.yearqtr(run.ho + hh/4 + 1e-8))
  fcst.out[[fcst.label]] <- fcst.smry
  fcst.draws[[fcst.label]] <- fcst
}

print(lapply(fcst.out, round, digits = 2))

# ---------------------------------------------------------------------------------------------
# released values of the targets over the past 4 quarters and the predictive densities
hist_y <- window(Yraw, start = tsp(Yraw)[2] - 3/4)
hist_df <- data.frame("date" = rep(as.numeric(time(hist_y)), N),
                      "target" = factor(rep(sl.target, each = nrow(hist_y)), levels = sl.target),
                      "value" = c(hist_y))

# draws of the predictive distributions at their target quarters (names of fcst.draws are ordered as grid.hh)
fcst_df <- do.call(rbind, lapply(names(fcst.draws), function(ll){
  data.frame("type" = ll, "target" = rep(sl.target, each = nrow(fcst.draws[[ll]])), "value" = c(fcst.draws[[ll]]))
}))
fcst_df$type <- factor(fcst_df$type, levels = names(fcst.draws))
fcst_df$target <- factor(fcst_df$target, levels = sl.target)
fcst_df$date <- (run.ho + grid.hh/4)[as.integer(fcst_df$type)]
fcst_df <- fcst_df[!is.na(fcst_df$value),]

# draws within the 1-99 percentiles to avoid long tails of the violins
viol_df <- do.call(rbind, lapply(split(fcst_df, list(fcst_df$type, fcst_df$target), drop = TRUE), function(dd){
  qq <- quantile(dd$value, c(0.01, 0.99))
  dd[dd$value >= qq[1] & dd$value <= qq[2],]
}))

# medians, 68 and 90 percent credible sets and top of the violins (for the labels)
med_df <- do.call(rbind, lapply(split(fcst_df, list(fcst_df$type, fcst_df$target), drop = TRUE), function(dd){
  qq <- quantile(dd$value, c(0.01, 0.05, 0.16, 0.5, 0.84, 0.95, 0.99))
  data.frame("type" = dd$type[1], "target" = dd$target[1], "date" = dd$date[1], "median" = qq[[4]],
             "q05" = qq[[2]], "q16" = qq[[3]], "q84" = qq[[5]], "q95" = qq[[6]], "top" = qq[[7]])
}))

# line from the last released value through the medians
path_df <- rbind(hist_df[hist_df$date == max(hist_df$date),], med_df[, c("date","target","median")] |> setNames(c("date","target","value")))

date_grid <- sort(unique(c(hist_df$date, med_df$date)))

pp_hist <- ggplot(hist_df, aes(x = date, y = value)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_violin(data = viol_df, aes(x = date, y = value, group = date), width = 0.18, fill = "grey85", colour = "grey30", linewidth = 0.25) +
  geom_linerange(data = med_df, aes(x = date, ymin = q05, ymax = q95, linewidth = "90% credible set"), inherit.aes = FALSE, colour = "grey30") +
  geom_linerange(data = med_df, aes(x = date, ymin = q16, ymax = q84, linewidth = "68% credible set"), inherit.aes = FALSE, colour = "grey30") +
  geom_line(data = path_df, linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  geom_point(data = med_df, aes(x = date, y = median, shape = "Median"), fill = "white", size = 1.6) +
  geom_text(data = med_df, aes(x = date, y = top, label = sprintf("%.2f", median)), vjust = -0.6, size = 3) +

  scale_linewidth_manual(values = c("68% credible set" = 1.6, "90% credible set" = 0.4), name = NULL) +
  scale_shape_manual(values = c("Median" = 21), name = NULL) +
  scale_x_continuous(breaks = date_grid, labels = format(as.yearqtr(date_grid + 1e-8), "%Y Q%q")) +
  scale_y_continuous(breaks = scales::breaks_extended(n = 10), expand = expansion(mult = c(0.05, 0.12))) +
  facet_wrap(. ~ target, ncol = 1, scales = "free_y") +

  xlab("") + ylab("Annualized growth rate (%)") +

  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.25),
        panel.spacing = unit(8, "pt"),
        axis.ticks = element_line(colour = "grey30", linewidth = 0.25),
        axis.ticks.length = unit(-2.5, "pt"),
        axis.text.x = element_text(margin = margin(t = 4)),
        axis.text.y = element_text(margin = margin(r = 4)),
        strip.text = element_text(colour = "black"))
print(pp_hist)
