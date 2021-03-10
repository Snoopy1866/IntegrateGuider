IntegrateGuider <- function()
{
    #-------------------------------------------------------
    #将输入的字符串转化为Numeric
    InputAs.Numeric <- function(chr)
    {
        n <- length(chr)
        value <- rep(NA,n)
        for(i in 1:n)
        {
            value[i] <- eval(parse(text = chr[i]))
        }
        
        value
    }
    #-------------------------------------------------------
    #将输入的字符串转化为Numeric Vector,支持3种输入模式：c(1,2,3) or (1,2,3) or 1,2,3
    InputAs.Vector <- function(chr)
    {
        
        n <- nchar(chr)
        if (substr(chr,1,2) == "c(")
        {
            InputAs.Numeric(strsplit(substr(chr,3,n - 1),",")[[1]])
            
        }
        else if (substr(chr,1,1) == "(")
        {
            InputAs.Numeric(strsplit(substr(chr,2,n - 1),",")[[1]])
        }
        else
        {
            InputAs.Numeric(strsplit(chr,",")[[1]])
        }
    }
    #-------------------------------------------------------
    #将输入的字符串转化为Function类型
    InputAs.Function <- function(chr)
    {
        if(exists(chr) && is.function(get(chr)))
        {
            f <- get(chr)
        }
        else
        {
            f <- function(x) 
            {
                expr <- parse(text = chr)
                Value <- function(x) eval(expr)
                Value(x)
            }
        }
        f
    }
    #-------------------------------------------------------
    InputPrompt <- function()
    {
        cat("请输入被积函数:( 例如:x^2*exp(-x^2) )\n")
        f <- InputAs.Function(readline())
        cat("请输入积分下限:\n")
        a <- InputAs.Numeric(readline())
        cat("请输入积分上限:\n")
        b <- InputAs.Numeric(readline())
        return(list(f = f,a = a,b = b))
    }
    #-------------------------------------------------------
    #判断积分结果
    OutputResult <- function(result)
    {
        if (result < -10e13)
            stop("积分结果可能为负无穷大")
        else if (result > 10e13)
            stop("积分结果可能为正无穷大")
        else if (abs(result) < 10e-13)
            return(0)
        else
            return(result)
    }
    #-------------------------------------------------------
    cat("欢迎使用数值积分向导!\n请选择积分类型:\n1. 定积分\n2. 重积分\n")
    integrate.type <- readline()
    if (integrate.type == "1")
    {
        cat("请选择积分方法(如果不输入将默认使用高斯求积公式):\n1. 高斯求积公式(可计算无穷区间的定积分)\n2. 龙贝格算法\n3. 复合辛普森\n4. 复合梯形\n5. 牛顿-科特斯\n")
        definite.integerate.method <- readline()
        if (any(definite.integerate.method == c("1","")))
        {
            cat("请输入被积函数(例如：x^2*exp(-x^2):\n")
            f <- InputAs.Function(readline())
            cat("请输入积分下限(-Inf 代表负无穷大,Inf 代表正无穷大):\n")
            a <- InputAs.Numeric(readline())
            cat("请输入积分上限(-Inf 代表负无穷大,Inf 代表正无穷大):\n")
            b <- InputAs.Numeric(readline())
            cat("需要使用几个高斯点?(默认为8)\n")
            chr <- readline()
            n <- ifelse(chr == "",8,as.numeric(chr))
            if (!(is.infinite(a) || is.infinite(b)))
            {
                cat("请选择具体的高斯求积公式(默认为勒让德求积公式):\n1. 勒让德求积公式\n2. 切比雪夫求积公式\n")
                gauss.method <- readline()
                if (any(gauss.method == c("1","")))
                {
                    result <- GaussIntegrate(f,a,b,n,gauss.method = "Legendre")
                }
                else if (gauss.method == "2")
                {
                    result <- GaussIntegrate(f,a,b,n,gauss.method = "Chebyshev")
                }
                else
                {
                    stop("输入有误!")
                }
            }
            else if (any(a == c(-Inf,Inf)) || any(b == c(-Inf,Inf)))
            {
                result <- GaussIntegrate(f,a,b,n)
            }
        }
        else if (definite.integerate.method == "2")
        {
            list <- InputPrompt()
            f <- list$f
            a <- list$a
            b <- list$b
            if (is.infinite(a) || is.infinite(b))
            {
                cat("无法使用龙贝格求积算法对无穷区间进行积分!\n改用高斯求积算法......\n")
                result <- GaussIntegrate(f,a,b,8)
            }
            else
            {
                cat("是否指定积分结果精度?请直接输入所需精度(默认为1e-5)\n")
                chr <- readline()
                eps <- ifelse(chr == "",1e-5,InputAs.Numeric(chr))
                cat("是否指定最大迭代次数?请直接输入最大迭代次数(默认为50次)\n注意:迭代次数过高可能导致运行时间过长...\n")
                chr <- readline()
                it_max <- ifelse(chr == "",50,InputAs.Numeric(chr))
                result <- Romberg(f,a,b,eps,it_max)
            }
        }
        else if (definite.integerate.method == "3")
        {
            cat("请输入被积函数:( 例如:x^2*exp(-x^2) )\n")
            f <- readline()
            expr <- parse(text = f)
            cat("请输入积分下限:\n")
            a <- InputAs.Numeric(readline())
            cat("请输入积分上限:\n")
            b <- InputAs.Numeric(readline())
            if (is.infinite(a) || is.infinite(b))
            {
                cat("无法使用复合辛普森求积算法对无穷区间进行积分!\n改用高斯求积算法......\n")
                result <- GaussIntegrate(InputAs.Function(f),a,b,8)
            }
            else
            {
                cat("请指定一个精度(默认为1e-5):\n")
                chr <- readline()
                eps <- ifelse(chr == "",1e-5,InputAs.Numeric(chr))
                result <- Simpson(expr,InputAs.Function(f),a,b,eps)
                return(cat("区间个数",result$区间个数,"\n积分结果",OutputResult(result$积分结果),"\n"))
            }
            
        }
        else if (definite.integerate.method == "4")
        {
            cat("请输入被积函数:( 例如:x^2*exp(-x^2) )\n")
            f <- readline()
            expr <- parse(text = f)
            cat("请输入积分下限:\n")
            a <- InputAs.Numeric(readline())
            cat("请输入积分上限:\n")
            b <- InputAs.Numeric(readline())
            if (is.infinite(a) || is.infinite(b))
            {
                cat("无法使用复合梯形求积算法对无穷区间进行积分!\n改用高斯求积算法......\n")
                result <- GaussIntegrate(InputAs.Function(f),a,b,8)
            }
            else
            {
                cat("请指定一个精度(默认为1e-5):\n")
                chr <- readline()
                eps <- ifelse(chr == "",1e-5,InputAs.Numeric(chr))
                result <- Trapezoid(expr,InputAs.Function(f),a,b,eps)
                return(cat("区间个数",result$区间个数,"\n积分结果",OutputResult(result$积分结果),"\n"))
            }
            
        }
        else if (definite.integerate.method == "5")
        {
            list <- InputPrompt()
            f <- list$f
            a <- list$a
            b <- list$b
            if (is.infinite(a) || is.infinite(b))
            {
                cat("无法使用牛顿-科特斯求积算法对无穷区间进行积分!\n改用高斯求积算法......\n")
                result <- GaussIntegrate(f,a,b,8)
            }
            else
            {
                cat("请输入要使用第几阶科特斯系数?(默认为4)\n")
                chr <- readline()
                n <- ifelse(chr == "",4,as.numeric(chr))
                if (n > 7 || n < 1)
                {
                    cat("指定的阶数不在1到7之间!\n将改用高斯求积算法......\n")
                    result <- GaussIntegrate(f,a,b,8)
                }
                else
                {
                    result <- NewtonCotes(f,a,b,n)
                }
            }
        }
        else
        {
            stop("输入有误!")
        }
    }
    else if (integrate.type == "2")
    {
        cat("请输入被积函数(被积函数必须事先定义，之后只需输入f即可):\n")
        f <- get(readline())
        cat("请设置随机数向量个数(默认为1e5)\n注意: 随机向量个数过大可能导致运算时间过长!\n")
        n <- readline()
        n <- ifelse(n == "",1e5,InputAs.Numeric(n))
        cat("请输入在此积分区域下自变量所能达到的下限构成的向量(必须为常数向量):\n")
        a <- InputAs.Vector(readline())
        cat("请输入在此积分区域下自变量所能达到的上限构成的向量(必须为常数向量):\n")
        b <- InputAs.Vector(readline())
        cat("积分区域是否为矩形区域?(即是否上下限均为常数向量)\n1. 是\n2. 否\n")
        is.cubic <- readline()
        if (is.cubic == "1")
        {
            cat("是否使用拟蒙特卡洛方法(默认,需要指定一个基数向量)?\n1. 是\n2. 否\n")
            is.QMC <- readline()
            if (any(is.QMC == c("1","")))
            {
                cat("请指定一个基数向量(基数必须为质数):\n")
                base <- readline()
                if (base == "")
                {
                    cat("您未指定基数向量,将改用蒙特卡洛方法...\n")
                    result <- MCorQMC(f,n,a = a,b = b,base = NULL,method = "MC")
                }
                else
                {
                    base <- InputAs.Vector(base)
                    result <- MCorQMC(f = f,n = n,a = a,b = b,base = base,method = "QMC")
                }
            }
            else if (is.QMC == "2")
            {   cat("请设置一个种子数(默认为1000):\n")
                seed <- readline()
                seed <- ifelse(seed == "",1000,InputAs.Numeric(seed))
                result <-MCorQMC(f,n,set_seed = seed,a = a,b = b,base = NULL,method = "MC")
            }
            else
            {
                stop("输入有误")
            }
        }
        else if (is.cubic == "2")
        {
            cat("请输入在此积分区域下自变量下限的约束条件所构成的向量(自变量真实的上下限构成的向量):\n")
            fail1 <- readline()
            if (is.null(fail1))
                stop("未输入任何值!")
            else
                fail1 <- get(fail1)
            cat("请输入在此积分区域下自变量上限的约束条件所构成的向量(自变量真实的上下限构成的向量):\n")
            fail2 <- readline()
            if (is.null(fail2))
                stop("未输入任何值!")
            else
                fail2 <- get(fail2)
            cat("是否使用拟蒙特卡洛方法(默认,需要指定一个基数向量)?\n1. 是\n2. 否\n")
            is.QMC <- readline()
            if (any(is.QMC == c("1","")))
            {
                cat("请指定一个基数向量(基数必须为质数):\n")
                base <- readline()
                if (base == "")
                {
                    cat("您未指定基数向量,将改用蒙特卡洛方法...\n")
                    result <- MCorQMC(f,n,a = a,b = b,fail1 = fail1,fail2 = fail2,base = NULL,method = "MC")
                }
                else
                {
                    base <- InputAs.Vector(base)
                    result <- MCorQMC(f,n,a = a,b = b,fail1 = fail1,fail2 = fail2,base = base,method = "QMC")
                }
            }
            else if (is.QMC == "2")
            {   cat("请设置一个种子数(默认为1000):\n")
                seed <- readline()
                seed <- ifelse(seed == "",1000,InputAs.Numeric(seed))
                result <- MCorQMC(f,n,set_seed = seed,a,b,fail1 = fail1,fail2 = fail2,base = NULL)
            }
            else
            {
                stop("输入有误")
            }
        }
        else
        {
            stop("输入有误!\n")
        }
        
    }
    else
    {
        stop("输入有误!\n")
    }
    return(cat("积分结果:",OutputResult(result),"\n"))
}


source('../IntegrateandDeriv/GaussIntegrate.R')
#调用高斯求积公式

source('../IntegrateandDeriv/Romberg.R')
#调用龙贝格求积公式

source('../Guider/Trapezoid.R')
#调用复合辛普森和复合梯形求积公式

source('../IntegrateandDeriv/NewtonCotes.R')
#调用牛顿-科特斯求积公式

source('../IntegrateandDeriv/MCorQMC.R')
#调用蒙特卡洛与拟蒙特卡洛公式

