library(countreg)
library(data.table)
library(dplyr)
library(flexmix)
library(ggplot2)

# This assumes there is a directory called "Data" in the working directory
dat <- data.table::fread("Data/MutSummary.csv")

p <- ggplot(dat, aes(x=Score, 
                     y=..density.., 
                     group=as.factor(Lineage),
                     fill=as.factor(Lineage),
                     colour=as.factor(Lineage))) +
    geom_histogram(position="identity")

printTable <- function(dat, filename, digits=3, hline_pos=0) {
    sink(filename)
    cat("\\begin{center}", '\n')
    dat %>%
        xtable::xtable(digits=digits) %>%
        print(
              floating=FALSE,
              include.rownames=FALSE,
              latex.environments="center",
              hline.after=c(0, hline_pos),
              sanitize.text.function=function(x){x}
             )
    cat("\\end{center}", '\n')
    sink()
}

getComponentParameters <- function(model,
                                   distribution
                                  ) {
    if(distribution == "poisson") {
        params <- parameters(model)
    } else if(distribution == "negative binomial") {
        mu <- model$coef
        theta <- summary(model)$theta
        params <- list(mu=mu,
                       theta=theta)
    }
    return(params)
}

getComponentDensity <- function(x,
                                params,
                                component_prop,
                                distribution
                               ) {
    dens <- NA
    if(distribution == "poisson") {
        dens <- component_prop*dpois(x, exp(params))
    } else if(distribution == "negative binomial") {
        dens <- component_prop*dnbinom(x, mu=exp(params$mu),
                                       size=params$theta
                                      )
    }
    return(dens)
}

getComponentProportions <- function(model,
                                    distribution
                                   ) {
    if(distribution == "poisson") {
    proportions <- attr(model, "prior")
    } else if (distribution == "negative binomial") {
        proportions <- 1
    }
    return(proportions)
}

getFullDensity <- function(x,
                           model,
                           distribution
                          ) {
    props <- model %>% getComponentProportions(distribution=distribution)
    params <- model %>% getComponentParameters(distribution=distribution)
    if(distribution == "negative binomial") {
        params <- list(params)
    }
    dens <- mapply(function(y, z) {
                          getComponentDensity(x, y, z,
                                              distribution=distribution)
                      },
                      params,
                      props
                     )
    if(length(x) > 1) {
        dens <- dens %>% apply(., 1, sum)
    } else {
        dens <- dens %>% sum
    }
    
    return(dens)
}

getNegativeBinomialModel <- function(values,
                                     ...
                                    ) {
    model <- glm.nb(values ~ 1)
    return(model)
}

getPoissonModel <- function(values,
                            k
                           ) {
    model <- flexmix(values ~ 1,
                     data=dat,
                     k=k,
                     model=FLXMRglm(family="poisson")
                    )
    return(model)
} 

getModelAIC <- function(model,
                        distribution) {
    if(distribution == "poisson") {
        aic <- attr(summary(model), "AIC") 
    } else if (distribution == "negative binomial") {
        aic <- model$aic
    }
    return(aic)
}

getModelBIC <- function(model) {
    return(model %>% BIC)
}

getPoissonMixtureDispersion <- function(model) {
    p <- model %>% getComponentProportions(distribution="poisson")
    lambda <- model %>% 
        getComponentParameters(distribution="poisson") %>%
        exp

    mu_M <- sum(p*lambda)
    mu2_M <- sum(p*lambda*(lambda + 1))

    dispersion <- mu2_M/mu_M - mu_M
    return(dispersion)
}

getNegativeBinomialDispersion <- function(model) {
    params <- model %>% getComponentParameters(distribution="negative binomial")
    lambda <- params$mu %>% exp
    theta <- params$theta
    dispersion <- 1 + lambda/theta
    return(dispersion)
}


getDensityMaximum <- function(observations,
                              model,
                              distribution
                             ) {
    max <- 0
    for(x in range(observations)) {
        val <- getFullDensity(x,
                              model,
                              distribution
                             )
        if(val > max) {
            max <- val
        }
    }
    return(max)
}

plotPoissonModel <- function(observations,
                      model,
                      filename,
                      distribution
                     ) {
    k <- attr(model, "k")
    params <- model %>% parameters
    props <- model %>% getComponentProportions(distribution=distribution)

    pdf(filename,
        width=8,
        height=8)
    par(mar=c(5,6,4,1)+.1)
    hist(observations, 
         prob=TRUE, 
         breaks=20,
         cex.lab=2,
         cex.axis=2,
         cex.main=2,
         ylim=c(0, 0.05),
         xlab="Mutation score",
         ylab="Frequency",
         main="Histogram and Poisson mixture p.m.f., K = " %>% paste0(k)
    )

    for(i in 1:k) {
    curve(getComponentDensity(x,
                              params=ifelse(distribution == "poisson",
                                        params[i],
                                        params[, i]),
                              component_prop=props[i],
                              distribution="poisson"
                             ),
         from=min(observations),
         to=max(observations),
         n=max(observations) - min(observations) + 1,
         add=TRUE,
         col="blue",
         lty=2,
         lwd=2
        )
    }
    curve(getFullDensity(x,
                         model,
                         distribution="poisson"
                        ),
         from=min(observations),
         to=max(observations),
         n=max(observations) - min(observations) + 1,
         add=TRUE,
         col="red"
        )

    legend("topleft",
           legend=c("Component density",
             "Full mixture density"),
           col=c("blue", "red"),
           lty=c(2, 1),
           cex=1.5
          )
    dev.off()
}

plotNegativeBinomialModel <- function(observations,
                                      model,
                                      filename
                                     ) {
    k <- 1
    params <- model %>% 
        getComponentParameters(distribution="negative binomial")

    pdf(filename,
        width=8,
        height=8
       )
    par(mar=c(5,6,4,1)+.1)
    hist(observations,
         prob=TRUE,
         breaks=20,
         cex.lab=2,
         cex.axis=2,
         cex.main=2,
         ylim=c(0, 0.05),
         xlab="Mutation score",
         ylab="Frequency",
         main="Histogram and negative binomial p.m.f."
    )
    curve(getFullDensity(x,
                         model,
                         distribution="negative binomial"
                        ),
          from=min(observations),
          to=max(observations),
          n=max(observations) - min(observations) + 1,
          add=TRUE,
          col="red"
         )
    legend("topleft",
           legend="Negative binomial density",
           col=c("red"),
           cex=1.5,
           lty=1
          )
    dev.off()
}

runAnalysis <- function(observations,
                        K,
						analysis_name
                       ) {
    pois_models <- list()
    pois_AICs <- {}
    pois_BICs <- {}
    pois_dispersions <- {}
    for(k in 1:K) {
        model <- getPoissonModel(observations, k)
        if(attr(model, "k") == k) {
            pois_models[[k]] <- model
            pois_AICs[k] <- getModelAIC(model,
                                        distribution="poisson"
                                       )
            pois_BICs[k] <- model %>% getModelBIC
            pois_dispersions[k] <- model %>% getPoissonMixtureDispersion
                                        
        }
    }
    
    nb_model <- getNegativeBinomialModel(observations)
    nb_AIC <- getModelAIC(nb_model,
                          distribution="negative binomial"
                         )
    nb_BIC <- nb_model %>% getModelBIC
    nb_dispersion <- nb_model %>% getNegativeBinomialDispersion
     
	figure_dir <- file.path("Figures",
							analysis_name
						   )
    unlink(figure_dir, rec=T)
	dir.create(figure_dir, showWarnings=FALSE)
    for(model in pois_models) {
        plotPoissonModel(observations,
                         model,
                         filename=file.path(figure_dir,
                                            paste0("pois_mixture_", 
                                            attr(model, "k"), ".pdf")
                                           ),
                         distribution="poisson"
                        )
    }

    plotNegativeBinomialModel(observations,
                              nb_model,
                              filename=file.path(figure_dir,
               				                     paste0("nb_mixture.pdf")
               				                    )
               				 )
    summary_df <- data.table(Model=paste("Poission, k = ", 1:length(pois_models)), 
                         AIC=pois_AICs,
                         BIC=pois_BICs,
                         "$\\estim D$"=pois_dispersions
                        )
    summary_df <- rbind(summary_df, 
                    data.table(Model="Negative Binomial", 
                               AIC=nb_AIC,
                               BIC=nb_BIC,
                               "$\\estim D$"=nb_dispersion
                              )
                   )

	table_dir <- file.path("Tables", analysis_name)
	dir.create(table_dir, showWarnings=FALSE)
	printTable(summary_df, file.path(table_dir,
								 "summary.tex"
								)
		      )
    return(summary_df)
}

getDaughterSums <- function(dat) {
    d1_scores <- dat[dat$SegregationID == 1, ]$Score
    d2_scores <- dat[dat$SegregationID == 2, ]$Score
    d_scores <- d1_scores + d2_scores
    return(d_scores)
}

getMotherSums <- function(dat) {
    m1_scores <- dat[dat$SegregationID == 3, ]$Score
    m2_scores <- dat[dat$SegregationID == 4, ]$Score
    m_scores <- m1_scores + m2_scores
    return(m_scores)
}

getFullSums <- function(dat) {
    d_scores <- dat %>% getDaughterSums
    m_scores <- dat %>% getMotherSums
    full_scores <- d_scores + m_scores
    return(full_scores)
}

K <- 4
# Run analysis for Mother 2 (M2)
runAnalysis(dat[dat$SegregationID == 4, ]$Score,
            K=K,
			analysis_name="M2"
           )

# Run daughter-only analysis
daughter_sums <- dat %>% getDaughterSums
runAnalysis(daughter_sums[daughter_sums > 80],
            K=K,
			analysis_name="Daughter"
           )
runAnalysis(daughter_sums,
            K=K,
			analysis_name="Daughter_full"
           )

# Run mother-only analysis
mother_sums <- dat %>% getMotherSums
runAnalysis(mother_sums[mother_sums > 60],
            K=K,
			analysis_name="Mother"
           )
runAnalysis(mother_sums,
            K=K,
			analysis_name="Mother_full"
           )

# Run analysis for all data
full_sums <- dat %>% getFullSums
runAnalysis(full_sums[full_sums > 170],
			K=K,
			analysis_name="all"
           )
runAnalysis(full_sums,
			K=K,
			analysis_name="all_full"
           )

# Paired t-test
m_sums <- dat %>% getMotherSums
d_sums <- dat %>% getDaughterSums
pdf("Figures/m_d_box.pdf")
boxplot(m_sums,
        d_sums,
        names=c("Mothers", "Daughters"),
        xlab="Cell type",
        ylab="Score",
        main="Mutation scores for mothers vs daughters"
       )
dev.off()
t_test <- t.test(m_sums, d_sums,
                 paired=TRUE)
       
