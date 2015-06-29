library(nlme)
library(qtl)

cimld <- function(pheno.file, geno.file, method=c("cim","sma"), n.perm=0)
{
    method <- match.arg(method)
    cross <- cimld.read(pheno.file, geno.file)

    pheno.names <- names(cross$doe)[-(1:3)]
    if (names(cross$doe)[1] == "ENVIR")
        pheno.names <- names(cross$doe)[-(1:4)]

    lodscores <- thresholds <- results <- NULL

    for (name in pheno.names) {
        cat("# Response:", name, "\n")
        varcomp <- cimld.varcomp(cross$doe, name)

        if (method == "cim") {
            if (!("draws" %in% names(cross$geno[[1]])))
                cross <- sim.geno(cross, step=1)
            z <- cimld.cim(cross, varcomp)
        }
        else {
            z <- cimld.sma(cross, varcomp)
        }

        if (is.null(lodscores))
            lodscores <- z
        else
            lodscores <- cbind(lodscores, z[,-(1:2),drop=FALSE])

        threshold <- 2.5

        if (n.perm > 0) {
            cat("# Permutation Test ...\n")
            y <- cross$doe[[name]]
            ycopy <- y
            x <- matrix(nrow=n.perm, ncol=ncol(z)-2)
            colnames(x) <- colnames(z)[-c(1,2)]

            for (i in 1:n.perm) {
                o <- sample(length(y))
                y <- y[o]
                cross$doe[[name]] <- y
                if (method == "cim")
                    temp <- cimld.cim(cross, varcomp)
                else
                    temp <- cimld.sma(cross, varcomp)
                x[i,] <- apply(temp[,-c(1,2),drop=FALSE], 2, max)
                if (i %% 10 == 0)
                    cat("  ",i,"\n",sep="")
            }

            zp <- apply(x, 2, quantile, c(0.99,0.95,0.9), na.rm=TRUE)
            rownames(zp) <- c("%1", "5%", "10%")
            if (zp[2,ncol(zp)] > threshold)
                threshold <- zp[2,ncol(zp)]
            thresholds <- cbind(thresholds, zp)

            cross$doe[[name]] <- ycopy
        }

        zq <- cimld.qtl(z, threshold=threshold)
        if (is.null(zq))
            next

        qtl <- makeqtl(cross, zq$chr, zq$pos)
        fit <- cimld.fit(cross, varcomp, qtl)
        results <- rbind(results, cbind(zq, fit))
    }

    if (!is.null(lodscores))
        write.csv(lodscores, file="cimld_output_lodscore.csv", quote=FALSE)

    if (!is.null(thresholds))
        write.csv(thresholds, file="cimld_output_threshold.csv", quote=FALSE)

    if (!is.null(results))
        write.csv(results, file="cimld_output_result.csv", quote=FALSE)

    cat("# CIMLD finished successfully!\n")
}

cimld.read <- function(pheno.file, geno.file)
{
    pheno <- read.table(pheno.file, header=TRUE, sep=",")
    stopifnot(!anyNA(pheno))

    cn <- toupper(names(pheno))
    j <- if (cn[1] != "ENVIR") 1 else 2
    stopifnot(cn[j] == "GROUP")
    stopifnot(cn[j+1] == "BLOCK")
    stopifnot(cn[j+2] == "LINE")

    cross <- read.cross("csv", file=geno.file, crosstype="riself")
    stopifnot(ncol(cross$pheno) == 1)

    obs <- match(pheno[,j+2], cross$pheno[,1])
    stopifnot(length(obs) > 10)

    names(pheno)[1:(j+2)] <- cn[1:(j+2)]

    pheno <- within(pheno, {
        GROUP <- factor(GROUP)
        BLOCK <- factor(BLOCK)
        LINE <- factor(LINE)
        BLOCK <- factor(GROUP:BLOCK)
    })

    if (j != 1) {
        pheno <- within(pheno, {
            ENVIR <- factor(ENVIR)
            BLOCK <- factor(ENVIR:BLOCK)
        })
    }

    cross$doe <- pheno
    cross$obs <- obs

    cross
}

cimld.sma <- function(cross, varcomp)
{
    y <- varcomp$w %*% cross$doe[[varcomp$name]]
    ac <- varcomp$w %*% varcomp$ac
    rss0.full <- sum(lm.fit(ac,y)$residuals^2)
    cnst.full <- 0.5*length(y)

    ic <- varcomp$ic
    n.lod <- 1

    header <- varcomp$name
    if (!is.null(ic)) {
        n.lod <- 3
        header <- c(paste(varcomp$name,c("G","GxE"),sep="."), header)
    }

    n.mar <- totmar(cross)
    z <- matrix(nrow=n.mar, ncol=n.lod)

    g <- pull.geno(cross) - 1

    for (j in 1:n.mar) {
        x <- g[cross$obs, j]
        if (anyNA(x)) {
            ind <- !is.na(x)
            yv <- y[ind]
            x0 <- ac[ind,]
            rss0 <- sum(lm.fit(x0,yv)$residuals^2)
            cnst <- 0.5*length(yv)
            if (is.null(ic)) {
                x1 <- varcomp$w[ind,ind] %*% x[ind]
                rss1 <- sum(lm.fit(cbind(x0,x1),yv)$residuals^2)
                z[j] <- cnst * log10(rss0/rss1)
            }
            else {
                w <- varcomp$w[ind,ind]
                xv <- x[ind]
                iv <- ic[ind,]
                x1 <- w %*% xv
                x2 <- w %*% model.matrix(~-1+iv:xv)
                rss1 <- sum(lm.fit(cbind(x0,x1),yv)$residuals^2)
                rss2 <- sum(lm.fit(cbind(x0,x1,x2),yv)$residuals^2)
                z[j,1] <- cnst * log10(rss0/rss1)
                z[j,2] <- cnst * log10(rss1/rss2)
                z[j,3] <- cnst * log10(rss0/rss2)
            }
        }
        else {
            x1 <- varcomp$w %*% x
            rss1 <- sum(lm.fit(cbind(ac,x1),y)$residuals^2)
            z[j] <- cnst.full * log10(rss0.full/rss1)
            if (!is.null(ic)) {
                x2 <- varcomp$w %*% model.matrix(~-1+ic:x)
                rss2 <- sum(lm.fit(cbind(ac,x1,x2),y)$residuals^2)
                z[j,2] <- cnst.full * log10(rss1/rss2)
                z[j,3] <- cnst.full * log10(rss0.full/rss2)
            }
        }
    }

    z <- as.data.frame(z)
    rownames(z) <- markernames(cross)
    colnames(z) <- header

    chr <- pos <- NULL
    for (i in 1:nchr(cross)) {
        map <- cross$geno[[i]]$map
        chr <- c(chr, rep(chrnames(cross)[i], length(map)))
        pos <- c(pos, as.numeric(map))
    }

    cbind(chr=chr, pos=pos, z)
}

cimld.cim <- function(cross, varcomp, window=10)
{
    window <- window / 2

    forw <- cimld.forwsel(cross, varcomp)

    mc <- NULL
    if (!is.null(forw))
        mc <- varcomp$w %*% (forw$geno[cross$obs,] - 1)

    y <- varcomp$w %*% cross$doe[[varcomp$name]]
    x0 <- varcomp$w %*% varcomp$ac
    rss0 <- sum(lm.fit(cbind(x0,mc),y)$residuals^2)
    cnst <- 0.5*length(y)

    ic <- varcomp$ic
    n.lod <- 1
    header <- varcomp$name
    if (!is.null(ic)) {
        n.lod <- 3
        header <- c(paste(varcomp$name,c("G","GxE"),sep="."), header)
    }

    results <- NULL

    for (i in 1:nchr(cross)) {
        draws <-  cross$geno[[i]]$draws - 1
        n.pos <- ncol(draws)
        n.draws <- dim(draws)[3]
        map <- attr(draws, "map")
        marnam <- names(map)

        adj.draws <- array(dim = c(length(y), n.pos, n.draws))
        for (k in 1:n.draws)
            adj.draws[,,k] <- varcomp$w %*% draws[cross$obs,,k]

        z <- matrix(nrow=n.pos, ncol=n.lod)
        zd <- matrix(nrow=n.draws, ncol=n.lod)

        for (j in 1:n.pos) {
            for (k in 1:n.draws) {
                x1 <- adj.draws[,j,k]
                rss1 <- sum(lm.fit(cbind(x0,mc,x1),y)$residuals^2)
                zd[k,1] <- cnst * log10(rss0/rss1)
                if (!is.null(ic)) {
                    g <- draws[cross$obs,j,k]
                    x2 <- varcomp$w %*% model.matrix(~-1+ic:g)
                    rss2 <- sum(lm.fit(cbind(x0,mc,x1,x2),y)$residuals^2)
                    zd[k,2] <- cnst * log10(rss1/rss2)
                    zd[k,3] <- cnst * log10(rss0/rss2)
                }
            }
            for (k in 1:n.lod)
                z[j,k] <- wtaverage(zd[,k])
        }

        if (!is.null(forw) && (chrnames(cross)[i] %in% forw$chrpos[,1])) {
            wh <- which(forw$chrpos[,1] == chrnames(cross)[i])
            lower.pos <- forw$chrpos[wh,2] - window
            upper.pos <- forw$chrpos[wh,2] + window
            for (j in 1:n.pos) {
                excl <- map[j] >= lower.pos & map[j] <= upper.pos
                if (!any(excl))
                    next
                ac <- mc[,-wh[excl]]
                rss0.new <- sum(lm.fit(cbind(x0,ac),y)$residuals^2)
                for (k in 1:n.draws) {
                    x1 <- adj.draws[,j,k]
                    rss1 <- sum(lm.fit(cbind(x0,ac,x1),y)$residuals^2)
                    zd[k,1] <- cnst * log10(rss0.new/rss1)
                    if (!is.null(ic)) {
                        g <- draws[cross$obs,j,k]
                        x2 <- varcomp$w %*% model.matrix(~-1+ic:g)
                        rss2 <- sum(lm.fit(cbind(x0,ac,x1,x2),y)$residuals^2)
                        zd[k,2] <- cnst * log10(rss1/rss2)
                        zd[k,3] <- cnst * log10(rss0.new/rss2)
                    }
                }
                for (k in 1:n.lod)
                    z[j,k] <- wtaverage(zd[,k])
            }
        }

        o <- grep("^loc-*[0-9]+", marnam)
        if (length(o) > 0)
            marnam[o] <- paste("c", names(cross$geno)[i], ".", marnam[o], sep="")

        z <- as.data.frame(z)
        colnames(z) <- header
        z <- cbind(chr=rep(chrnames(cross)[i],n.pos), pos=as.numeric(map), z)
        rownames(z) <- marnam

        results <- rbind(results, z)
    }

    results
}

cimld.forwsel <- function(cross, varcomp, maxsize=10, alpha=0.05)
{
    g <- pull.geno(cross)
    if (anyNA(g))
        g <- pull.geno(fill.geno(cross))

    x <- varcomp$w %*% (g[cross$obs,] - 1)

    y <- varcomp$w %*% cross$doe[[varcomp$name]]

    x0 <- varcomp$w %*% varcomp$ac
    fit <- lm.fit(x0, y)
    rss0 <- sum(fit$residuals^2)
    df0 <- fit$df.residual

    n.mar <- totmar(cross)
    bonf <- alpha / n.mar
    ignore <- logical(n.mar)
    chosen <- NULL

    while (length(chosen) < maxsize) {
        cur.pval <- 1
        cur.idx <- cur.rss <- cur.df <- NULL

        for (i in 1:n.mar) {
            if (ignore[i])
                next

            fit <- lm.fit(cbind(x0,x[,i]), y)
            rss1 <- sum(fit$residuals^2)
            df1 <- fit$df.residual

            if (df1 > 0 && rss1 > 0 && df0 > df1 && rss0 > rss1) {
                fval <- (rss0-rss1)/(df0-df1) / (rss1/df1)
                pval <- pf(fval, df0-df1, df1, lower.tail=FALSE)
                if (pval < cur.pval) {
                    cur.pval <- pval
                    cur.idx <- i
                    cur.rss <- rss1
                    cur.df <- df1
                }
            }
        }

        if (cur.pval > bonf)
            break

        chosen <- c(chosen, cur.idx)
        ignore[cur.idx] = TRUE

        x0 <- cbind(x0, x[,cur.idx])
        rss0 <- cur.rss
        df0 <- cur.df
    }

    if (length(chosen) == 0)
        return(NULL)

    mar <- markernames(cross)[chosen]
    chrpos <- find.markerpos(cross, mar)

    list(mar=mar, chrpos=chrpos, geno=g[,mar,drop=FALSE])
}

cimld.varcomp <- function(pheno, name, GxE=TRUE)
{
    if (names(pheno)[1] == "ENVIR")
        return(cimld.varcomp.envir(pheno, name, GxE))

    fmla <- "~GROUP+LINE"
    fit <- lme(as.formula(paste(name,fmla)), random=~1|BLOCK, data=pheno)

    cat("# Covariance Parameter Estimates\n")
    print(VarCorr(fit))
    cat("# Tests of Fixed Effects\n")
    print(anova(fit))

    vc <- as.numeric(VarCorr(fit))
    var.block <- vc[1]
    var.error <- vc[2]

    Z <- model.matrix(~-1+BLOCK, data=pheno)
    n <- nrow(Z)
    V <- tcrossprod(Z)*var.block + diag(var.error,n)

    V <- V / var.error
    V <- backsolve(chol(V), diag(n))
    V <- t(V)

    addcov <- model.matrix(~GROUP, data=pheno)
    intcov <- NULL

    list(vc=vc, w=V, ac=addcov, ic=intcov, name=name)
}

cimld.varcomp.envir <- function(pheno, name, GxE)
{
    fmla <- if (GxE) "~ENVIR+ENVIR:GROUP+LINE+LINE:ENVIR" else "~ENVIR+ENVIR:GROUP+LINE"
    fit <- lme(as.formula(paste(name,fmla)), random=~1|BLOCK, data=pheno)

    cat("# Covariance Parameter Estimates\n")
    print(VarCorr(fit))
    cat("# Tests of Fixed Effects\n")
    print(anova(fit))

    vc <- as.numeric(VarCorr(fit))
    var.block <- vc[1]
    var.error <- vc[2]

    Z <- model.matrix(~-1+BLOCK, data=pheno)
    n <- nrow(Z)
    V <- tcrossprod(Z)*var.block + diag(var.error,n)

    V <- V / var.error
    V <- backsolve(chol(V), diag(n))
    V <- t(V)

    addcov <- model.matrix(~ENVIR+ENVIR:GROUP, data=pheno)
    intcov <- if (GxE) model.matrix(~ENVIR, data=pheno)[,-1,drop=FALSE] else NULL

    list(vc=vc, w=V, ac=addcov, ic=intcov, name=name)
}

cimld.qtl <- function(d, threshold=2.5, window=10)
{
    stopifnot(ncol(d) == 3 || ncol(d) == 5)

    window <- window / 2
    results <- NULL

    for (chr in unique(d$chr)) {
        z <- d[d$chr==chr,]
        pos <- z[,2]
        lod <- z[,ncol(z)]

        qidx <- which(lod >= threshold)

        qdel <- NULL
        for (i in qidx) {
            lower <- pos[i] - window
            upper <- pos[i] + window
            if (any(lod[(pos>lower)&(pos<upper)] > lod[i]))
                qdel <- c(qdel, i)
        }
        qidx <- setdiff(qidx, qdel)

        qdel <- NULL
        for (i in which(diff(lod[qidx])==0))
            qdel <- c(qdel, qidx[i])
        qidx <- setdiff(qidx, qdel)

        lsi <- usi <- NULL
        for (i in qidx) {
            wh <- which(lod < lod[i]-1.5)
            j <- wh[tail(which(wh<i),n=1)]
            k <- wh[head(which(wh>i),n=1)]
            if (length(j) == 0)
                j <- i
            if (length(k) == 0)
                k <- i
            lsi <- c(lsi, pos[j])
            usi <- c(usi, pos[k])
        }

        qdel <- NULL
        for (i in which(diff(lsi)==0)) {
            if (usi[i] == usi[i+1])
                qdel <- c(qdel, (if (lod[qidx[i]]<lod[qidx[i+1]]) i else i+1))
        }
        if (!is.null(qdel)) {
            qidx <- qidx[-qdel]
            lsi <- lsi[-qdel]
            usi <- usi[-qdel]
        }

        if (length(qidx) != 0)
            results <- rbind(results, cbind(z[qidx,], lsi=lsi, usi=usi))
    }

    if (is.null(results))
        return(NULL)

    if (ncol(d) == 3)
        names(results)[3] <- "lod"
    else
        names(results)[3:5] <- paste(c("G.","GxE.",""),"lod",sep="")

    rownames(results) <- paste(tail(names(d),n=1), rownames(results), sep=".")

    results
}

cimld.fit <- function(cross, varcomp, qtl)
{
    y <- varcomp$w %*% cross$doe[[varcomp$name]]
    x0 <- varcomp$ac

    xe <- varcomp$ic
    if (!is.null(xe))
        xe[apply(xe,1,sum)==0,] <- -1

    rss0 <- sum((y-mean(y))^2)
    cnst <- 0.5*length(y)

    n.qtl <- qtl$n.qtl
    n.draws <- dim(qtl$geno)[3]

    lod <- freq <- coeff <- NULL

    for (k in 1:n.draws) {
        x1 <- qtl$geno[cross$obs,,k] - 1
        x1[x1==0] <- -1
        freq <- rbind(freq, apply(x1,2,sum))
        x2 <- NULL
        if (!is.null(xe)) {
            for (i in 1:n.qtl) {
                for (j in 1:ncol(xe))
                    x2 <- cbind(x2, x1[,i]*xe[,j,drop=FALSE])
            }
        }
        x <- varcomp$w %*% cbind(x0,x1,x2)
        fit <- lm.fit(x, y)
        rss1 <- sum(fit$residuals^2)
        lod <- c(lod, cnst*log10(rss0/rss1))
        coeff <- rbind(coeff, fit$coefficients[-(1:ncol(x0))])
    }

    # Method modified from R/qtl package, http://www.rqtl.org/.
    # Broman KW, Wu H, Sen Ś, Churchill GA (2003) R/qtl: QTL mapping in experimental crosses. Bioinformatics 19:889-890.

    wts <- lod * log(10)
    tot.wt <- wts[1]
    for (e in wts[-1])
        tot.wt <- addlog(tot.wt, e)
    wts <- exp(wts - tot.wt)

    ord <- order(lod)
    m = floor(0.5*log(n.draws)/log(2))
    idx <- c(head(ord,n=m),tail(ord,n=m))

    wts[idx] <- 0
    wts <- wts / sum(wts)

    a <- apply(matrix(wts,nrow=n.draws,ncol=ncol(coeff))*coeff, 2, sum)
    p <- apply(freq[-idx,],2,mean) / length(y)
    vg <- (1-p^2) * a[1:n.qtl]^2
    h <- vg / (sum(vg) + varcomp$vc[2])

    if (is.null(xe))
        return(data.frame(est=a, hsq=h))

    est.GxE <- matrix(a[-(1:n.qtl)], nrow=n.qtl, byrow=TRUE)
    colnames(est.GxE) <- paste("GxE.", tail(names(a),n=ncol(xe)), sep="")

    data.frame(est=a[1:n.qtl], hsq=h, est.GxE=est.GxE)
}

# Util functions modified from R/qtl package, http://www.rqtl.org/.
# Broman KW, Wu H, Sen Ś, Churchill GA (2003) R/qtl: QTL mapping in experimental crosses. Bioinformatics 19:889-890.

wtaverage <- function(x)
{
    n <- length(x)
    m <- floor(0.5*log(n)/log(2))
    x <- sort(x)
    x <- x[(m+1):(n-m)]
    mean(x) + 0.5*log(10)*var(x)
}

addlog <- function(a, b)
{
    if (b > a+200)
        return(b)
    if (a > b+200)
        return(a);
    return(a + log1p(exp(b-a)))
}
