DerivGuider <- function()
{
    #-------------------------------------------------------
    #将输入的字符串转化为Numeric Vector
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
    InputAs.Function <- function(chr,x)
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
    #三次样条求导法向导
    CubicSplineDerivGuider <- function(x,y)
    {
        cat("请选择边界条件!\n1. 已知两端的一阶导数值\n2. 已知两端的二阶导数值\n3. 原函数为周期函数\n")
        spline.condition <- readline()
        if (spline.condition == "1")
        {
            b <- 1
            cat("请输入两端的一阶导数值向量\n")
            bz <- InputAs.Vector(readline())
            cat("请输入所要求导数值的点的横坐标值\n")
            x.value <- InputAs.Numeric(readline())
            CubicSplineDeriv(x,y,x.value,b,bz)
        }
        else if (spline.condition == "2")
        {
            b <- 2
            cat("请输入两端的二阶导数值向量\n")
            bz <- InputAs.Vector(readline())
            cat("请输入所要求导数值的点的横坐标值\n")
            x.value <- InputAs.Numeric(readline())
            CubicSplineDeriv(x,y,x.value,b,bz)
        }
        else if (spline.condition == "3")
        {
            if (y[1] != y[length(y)])
                stop("周期函数两端点函数值必须相等!\n")
            else
            {
                b <- 3
                cat("请输入所要求导数值的点的横坐标值\n")
                x.value <- InputAs.Numeric(readline())
                CubicSplineDeriv(x,y,x.value,b)
            }
        }
    }
    #-------------------------------------------------------
    cat("欢迎使用数值微分向导!\n请选择微分的方法:\n1. 理查森外推法\n2. 插值型求导法\n")
    deriv.method <- readline()
    if (deriv.method == "1")
    {
        cat("请输入要求导的函数( 例如:x^2*exp(-x) ):\n")
        f <- InputAs.Function(readline())
        cat("请输入要求导的x的值:\n")
        x <- InputAs.Numeric(readline())
        cat("请输入初始步长(默认为0.1):\n")
        h <- readline()
        h <- ifelse(h == "",0.1,InputAs.Numeric(h))
        Waitui(f,x,h)
    }
    else if (deriv.method == "2")
    {
        cat("请输入插值点对应的x向量:\n")
        x <- InputAs.Vector(readline())
        cat("请输入插值点对应的y向量:\n")
        y <- InputAs.Vector(readline())
        if (length(x) == 2)
        {
            TwoPoint(x,y)
        }
        else if (length(x) == 3)
        {
            if (x[2] != mean(x))
            {
                cat("节点不是等距的,将使用三次样条求导法......\n")
                CubicSplineDerivGuider(x,y)
            }
            else
            {
                cat("您输入了三个点,是否使用三次样条求导法(默认)?\n1. 是\n2. 否\n")
                is.spline <- readline()
                if (any(is.spline == c("N","n","2")))
                {
                  ThreePoint(x,y)
                }
                if (any(is.spline == c("Y","y","1","")))
                {
                  CubicSplineDerivGuider(x,y)
                }
            }
        }
        else
        {
            cat("您输入了三个以上的点，将使用三次样条求导法!\n")
            CubicSplineDerivGuider(x,y)
        }
    }
}

source('../IntegrateandDeriv/MyDeriv.R')
#调用数值微分函数