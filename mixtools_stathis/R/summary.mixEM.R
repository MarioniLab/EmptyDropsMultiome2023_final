summary.mixEM <- function(object, digits=6, ...){
  x <- object
  o <- switch(x$ft,
              "multmixEM" = rbind(x$lambda, t(x$theta)),
              "normalmixEM" = rbind(x$lambda, x$mu, x$sigma),
              "repnormmixEM" = rbind(x$lambda, x$mu, x$sigma),
              "regmixEM" = rbind(x$lambda, x$sigma, x$beta),
              "regmixEM.lambda" = rbind(x$lambda, x$sigma, x$beta),
              "regmixEM.mixed" = rbind(x$lambda, x$sigma, x$beta),
              "regmixEM.loc" = rbind(x$sigma, x$beta),
              "regmixEM.chgpt" = rbind(x$lambda, x$sigma),
              "logisregmixEM" = rbind(x$lambda, x$beta),
              "poisregmixEM" = rbind(x$lambda, x$beta),
              "expRMM_EM" = rbind(x$lambda, x$rate),       # Reliability Mixt Model censored exp mixture
              "weibullRMM_SEM" = rbind(x$lambda, x$shape, x$scale), # RMM censored weibull mixture
              "mvnormalmixEM" = rbind(x$lambda, matrix(unlist(x$mu), byrow=FALSE, ncol=length(x$lambda))),
              "normalmixMMlc" = rbind(x$lambda, x$mu, x$sigma),
              stop("Unknown mixEM object of type ", x$ft))
  colnames(o) <- paste("comp",1:ncol(o))
  rownames(o) <- switch(x$ft,
                        "multmixEM" = c("lambda", paste("theta", 1:ncol(x$theta), sep="")),
                        "normalmixEM" = c("lambda", "mu", "sigma"),
                        "repnormmixEM" = c("lambda", "mu", "sigma"),
                        "regmixEM" = c("lambda", "sigma", paste("beta", 1:nrow(x$beta), sep="")),
                        "regmixEM.lambda" = c("lambda", "sigma", paste("beta", 1:nrow(x$beta), sep="")),
                        "regmixEM.mixed" = c("lambda", "sigma", paste("beta", 1:nrow(x$beta), sep="")),
                        "regmixEM.loc" = c("sigma", paste("beta", 1:nrow(x$beta), sep="")),
                        "regmixEM.chgpt" = c("lambda", "sigma"),
                        "logisregmixEM" = c("lambda", paste("beta", 1:nrow(x$beta), sep="")),
                        "poisregmixEM" = c("lambda", paste("beta", 1:nrow(x$beta), sep="")),
                        "expRMM_EM" = c("lambda", "rate"),
                        "weibullRMM_SEM" = c("lambda", "shape", "scale"),
                        "mvnormalmixEM" = c("lambda", paste("mu", 1:length(x$mu[[1]]), sep="")),
                        "normalmixMMlc" = c("lambda", "mu", "sigma"))
	cat("summary of", x$ft, "object:\n")
	print(o, digits=digits)
	cat("loglik at estimate: ", x$loglik, "\n")
	if (x$ft == "expRMM_EM") cat(100*(1-mean(x$d)), "% of the data right censored\n")
	if (x$ft == "weibullRMM_SEM") cat(100*(1-mean(x$d)), "% of the data right censored\n")
}

