this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
# setwd("E:/Users/17531/Documents/R/R.7/Guider")
source(file = "DerivGuider.R")
source(file = "IntegrateGuider.R")
options(digits = 10)
Guider <- function()
{
    cat("您想计算什么?\n1. 数值积分\n2. 数值微分\n")
    cat("注意:

        如果您需要计算重积分,请提前定义好被积函数与积分区域的约束条件
        例如:
            使用(拟)蒙特卡洛方法计算函数 f(x,y,z) = z^3*sin(y) 在
            区域{(x,y,z)|0 < x < 2,0 < y < pi,0 < z < cos(y)}上的积分。
            需先定义f <- function(x) x[3]^3*sin(x[2])
            (其中变量x,y,z分别用x[1],x[2],x[3]表示)
            以及约束条件:
            fail1 <- function(x) c(0,0,0)
            fail2 <- function(x) c(2*pi,pi/2,cos(x[2]))
            对于积分区域上下限为常数的,您无需定义fail1与fail2


        如果您已经定义了所需函数,请直接输入对应的数字 1 或 2
        如果您未准备好所需函数，请按Enter结束,定义完所需函数后再执行此向导!\n")
    num <- readline()
    if (num == "1")
        IntegrateGuider()
    else if(num == "2")
        DerivGuider()
    else
        return(cat(""))
}
Guider()