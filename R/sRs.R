sRs <- function(x){

		#x <- Xdat[[5]]

			Nn <- dim(x)[1]
			Pp <- dim(x)[2]

			standX <- scale((x))
			useSVD <- svd(standX)
			#range(sapply(1:Nn, function(t) sum(useSVD$u[t,]^2*useSVD$d^4))/(Pp-1)^2)
			#sum(sapply(1:Nn, function(t) sum(useSVD$u[t,]^2)))

			adj4ListR <- sapply(1:Pp, function(t) sum(useSVD$v[t,-Nn]^2*useSVD$d[-Nn]^4))
			#adj4ListL <- sapply(1:Nn, function(t) sum(useSVD$u[t,]^2*useSVD$d^4))
			sqR <- adj4ListR/(Nn-1)^2; range(sqR);
			maxR <- useSVD$d[1]^2/(Nn-1)

			signal <- (maxR/Pp)
			#c(max(Pp, Nn-1)/(Nn-1), maxR)

			condition_num = useSVD$d[1]/useSVD$d[min(Pp, Nn-1)]
			singular_val <- useSVD$d
			eigenval <- useSVD$d^2/sum(useSVD$d^2)*min(Pp, Nn)
			RED = sqrt(sum(sqR)-Pp)/sqrt(Pp*(Pp-1)); RED^2

			singular_val_use <- singular_val[1:min(Pp, Nn-1)]



# ### local; strong or weak:

## use this to find outliers when p < n

if (Pp < Nn){


	 # # minRj <- sapply(1:Pp, function(t) sum(abs(useSVD$v[t,]))); range(minRj)
	 # aa <- c(sum(eigenval[eigenval>1])/eigenval[Pp-1], eigenval[1]/min(eigenval))
	 # aa/sum(aa)

	 # pickOutlier <- (minRj > (median(minRj) + 3*sd(minRj)))| (minRj < (median(minRj) - 3*sd(minRj)))

# # 	 # ww <- aa*c(sum(!pickOutlier), sum(pickOutlier))/sum(aa*c(sum(!pickOutlier), sum(pickOutlier)))
	 # # ww


	 # # sqRsL <- ww[2]*ifelse(sum(pickOutlier) > 0, sum(sqR[pickOutlier]-1)/sum(pickOutlier)/(eigenval[1]-1), 0); sqRsL
	 # # sqRsB <- ifelse(sum(!pickOutlier) > 0, sum(sqR[!pickOutlier]-1)/sum(!pickOutlier)/(Pp-1), 0); sqRsB

	 # # sqRs <- sqRsL + sqRsB; sqRs

	 w1 <- sum((min(Pp, Nn)/singular_val_use[singular_val_use > abs(sqrt(Pp)-sqrt(Nn))]^2))/sum((min(Pp, Nn)/singular_val_use^2))

     sqRsL <- (sum(sqR)-Pp)/(Pp*(singular_val[1]^2/(Nn-1)-1))*(1-w1); sqRsL
	 sqRsB <- (sum(sqR)-Pp)/(Pp*(Pp-1))*w1; sqRsB

	 sqRs <- sqRsL + sqRsB; sqRs

	 sqRsL <- (sum(sqR)-Pp)/(Pp*(singular_val[1]^2/(Nn-1)-1))*((1-w1) + sum(singular_val_use[singular_val_use^2 < Pp]^2)/sum(singular_val_use^2))/2; sqRsL
	 sqRsB <- (sum(sqR)-Pp)/(Pp*(Pp-1))*(w1+sum(singular_val_use[singular_val_use^2 > Pp]^2)/sum(singular_val_use^2))/2; sqRsB
	 sqRs <- sqRsL + sqRsB; sqRs

 } else {

# ## when p > n; outliers become indistinguishable

	 w1 <- sum((min(Pp, Nn)/singular_val_use[singular_val_use > abs(sqrt(Pp)-sqrt(Nn))]^2))/sum((min(Pp, Nn)/singular_val_use^2))

	 sqRsL <- (sum(sqR)-Pp)/(Pp*(singular_val[1]^2/(Nn-1)-1))*((1-w1) + sum(singular_val_use[singular_val_use^2 < Pp]^2)/sum(singular_val_use^2))/2; sqRsL
	 sqRsB <- (sum(sqR)-Pp)/(Pp*(Pp-1))*(w1+sum(singular_val_use[singular_val_use^2 > Pp]^2)/sum(singular_val_use^2))/2; sqRsB
	 sqRs <- sqRsL + sqRsB; sqRs

 }


	#sqRs_upper <- (sum(sqR)-Pp)/(Pp*(Pp-1));
	#sum(dist1 > dist2)/Pp*sqrt(sum(sqR[dist1 > dist2])/sum(dist1 > dist2)/Pp) + ifelse(sum(dist1 < dist2) > 0,  sum(dist1 < dist2)/Pp*sqrt(sum(sqR[dist1 < dist2])/sum(dist1 < dist2)/maxR), 0); sqRs
	# ifelse(sum(dist1 < dist2) > 0,  sum(dist1 < dist2)/Pp*(sum(sqR[dist1 < dist2]) - sum(dist1 < dist2))/sum(dist1 < dist2)/(sum(dist1 < dist2)-1), 0); sqRs

    result_list <- list(sqR, condition_num, eigenval, RED, singular_val, sqRs, sqRsL, sqRsB)

	#sqrt(c(eigenval[1]/median(eigenval), median(eigenval)/min(eigenval)))

	if (Nn > Pp) {
	eigen_vals <- ((useSVD$d)^2/Nn);
	#sum(eigen_vals^2)/sqrt(Nn*Pp^2)
	#sum((eigen_vals-1)^2)/sqrt(Nn*(Pp-1)^2)
	#sum((eigen_vals)^2); sum(useSVD$d^4)/Nn^2
	#sqrt(sum((eigen_vals-1)^2))/(Pp/sqrt(Pp-1))

	#sapply(1:Pp, function(t) sum(useSVD$v[t,]^2*useSVD$d^2))

	y <- stats::rnorm(Nn)
	datFrame <- as.data.frame(cbind(y, standX))
	names(datFrame)[1] <- "y"
	model1 <- stats::lm(y ~ ., data = datFrame)
	vif <- car::vif(model1)

	#condition_num <- sapply(1:Pp, function(jj){
	#	svd_stand <- svd(standX[,-jj])
	#	svd_stand$d[1]/svd_stand$d[Pp-1]
	#})
	result_list <- list(cbind(sqR, vif), condition_num, eigenval, RED, singular_val, sqRs, sqRsL, sqRsB)
	}

	return(result_list)
	}

