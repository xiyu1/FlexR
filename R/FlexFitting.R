#' Easy curve fitting in R
#' @import brooms
#' @import dplyr
#' @importFrom brooms tidy glance argument
#' @importFrom dplyr %>% mutate left_join select everything mutate_if
#' @description \itemize{ \item flex_fit(): Fast fitting of multiple equations. \item #'     show_formulalist(): Show the original equation
#' @param y A numerical vector
#' @param x A numerical vector
#' @param logistic a logical value indicating whether logistic formula should be
#' fitted.
#' @param aug a logical value indicating whether augment should be export
#' @return flex_fit(): Returns a tibble or list
#' @examples
#' # Compute a formula
#' flex_fit(mtcars$mpg,mtcars$disp)
#' @export

# function body
flex_fit <- function(y, x, logistic = FALSE, aug = FALSE){
  models <- list()#p值等
  if(all(is.numeric(x),is.numeric(y))){
    if(!logistic){
      if(all(x>0,y>0)){
        modlist <- list(Linear = lm(y~x),#线性
                        Quadratic = lm(y~x+I(x^2)),#二次
                        Compound = lm(log(y)~x),#复合
                        Growth = lm(log(y)~x),#生长
                        Logarithmic = lm(y~log(x)),#对数
                        S = lm(log(y)~I(1/x)),#S曲线
                        Exponential = lm(log(y)~x),#指数
                        Inverse = lm(y~I(1/x)),#反函数
                        Power = lm(log(y)~log(x))#幂函数
        )
        if(aug){
          models <- modlist %>% map(broom::augment)
        }else{
          models[['tidy']] <- modlist %>% map(broom::tidy)
          models[['glance']] <- modlist %>% map(broom::glance)
          #对系数进行还原
          models$tidy$Compound[,2] <- exp(models$tidy$Compound[,2])#复合
          models$tidy$Exponential[1,2] <- exp(models$tidy$Exponential[1,2])#指数
          models$tidy$Power[1,2] <- exp(models$tidy$Power[1,2])#幂函数
          #生成表达式过程-----
          #1系数小数点保留,为提取方程的系数做准备
          models[["tidy"]] <- models[["tidy"]] %>% map(function(df){
            round(df %>% dplyr::select_if(is.numeric),2)
          })
          #提取模型为数据框
          models$glance %>% bind_rows(.id = "model") -> mod
          #生成方程
          df <- data.frame(
            mod =  c(
              "Linear",
              "Quadratic",
              "Compound",
              "Growth",
              "Logarithmic",
              "S",
              "Exponential",
              "Inverse",
              "Power"
            ),
            formula = c(
              #线性
              ifelse(
                models[["tidy"]][["Linear"]][["estimate"]][2] >= 0,
                paste0('y=', models[["tidy"]][["Linear"]][["estimate"]][1], '+',
                       models[["tidy"]][["Linear"]][["estimate"]][2], 'x'),
                paste0('y=', models[["tidy"]][["Linear"]][["estimate"]][1] ,
                       models[["tidy"]][["Linear"]][["estimate"]][2], 'x')
              ),
              #二次
              if (all(models[["tidy"]][["Quadratic"]][["estimate"]][2:3] >=
                      0)) {
                paste0('y=', models[["tidy"]][["Quadratic"]][["estimate"]][1],
                       '+', models[["tidy"]][["Quadratic"]][["estimate"]][2], 'x',
                       '+', models[["tidy"]][["Quadratic"]][["estimate"]][3], 'x^2')
              } else if (models[["tidy"]][["Quadratic"]][["estimate"]][2] >=
                         0 &
                         models[["tidy"]][["Quadratic"]][["estimate"]][3] <=
                         0) {
                paste0('y=', models[["tidy"]][["Quadratic"]][["estimate"]][1], '+',
                       models[["tidy"]][["Quadratic"]][["estimate"]][2], 'x',
                       models[["tidy"]][["Quadratic"]][["estimate"]][3], 'x^2')
              } else if (models[["tidy"]][["Quadratic"]][["estimate"]][2] <=
                         0 &
                         models[["tidy"]][["Quadratic"]][["estimate"]][3] >=
                         0) {
                paste0('y=', models[["tidy"]][["Quadratic"]][["estimate"]][1],
                       models[["tidy"]][["Quadratic"]][["estimate"]][2], 'x', '+',
                       models[["tidy"]][["Quadratic"]][["estimate"]][3], 'x^2')
              } else if (models[["tidy"]][["Quadratic"]][["estimate"]][2] <=
                         0 &
                         models[["tidy"]][["Quadratic"]][["estimate"]][3] <= 0) {
                paste0('y=', models[["tidy"]][["Quadratic"]][["estimate"]][1],
                       models[["tidy"]][["Quadratic"]][["estimate"]][2], 'x',
                       models[["tidy"]][["Quadratic"]][["estimate"]][3], 'x^2')
              },
              #复合
              paste0('y=', models[["tidy"]][["Compound"]][["estimate"]][1] , '(',
                     models[["tidy"]][["Compound"]][["estimate"]][2], '^x', ')'),
              #生长
              ifelse(
                models[["tidy"]][["Growth"]][["estimate"]][2] >= 0,
                paste0('y=', 'exp(', models[["tidy"]][["Growth"]][["estimate"]][1] ,
                       '+', models[["tidy"]][["Growth"]][["estimate"]][2], 'x', ')'),
                paste0('y=', 'exp(', models[["tidy"]][["Growth"]][["estimate"]][1] ,
                       models[["tidy"]][["Growth"]][["estimate"]][2], 'x', ')')
              ),
              #对数
              ifelse(
                models[["tidy"]][["Logarithmic"]][["estimate"]][2] >= 0,
                paste0('y=', models[["tidy"]][["Logarithmic"]][["estimate"]][1] , '+',
                       models[["tidy"]][["Logarithmic"]][["estimate"]][2],
                       'ln(x)'),
                paste0('y=', models[["tidy"]][["Logarithmic"]][["estimate"]][1] ,
                       '', models[["tidy"]][["Logarithmic"]][["estimate"]][2],
                       'ln(x)')
              ),
              #S
              ifelse(
                models[["tidy"]][["S"]][["estimate"]][2] >= 0,
                paste0('y=', 'exp(', models[["tidy"]][["S"]][["estimate"]][1] , '+',
                       models[["tidy"]][["S"]][["estimate"]][2], '/x', ')'),
                paste0('y=', 'exp(', models[["tidy"]][["S"]][["estimate"]][1] ,
                       models[["tidy"]][["S"]][["estimate"]][2], '/x', ')')
              ),
              #指数函数
              ifelse(
                models[["tidy"]][["Exponential"]][["estimate"]][2] >= 0,
                paste0('y=', models[["tidy"]][["Exponential"]][["estimate"]][1] ,
                       'exp(', models[["tidy"]][["Exponential"]][["estimate"]][2], 'x',
                       ')'),
                paste0('y=', models[["tidy"]][["Exponential"]][["estimate"]][1] ,
                       'exp(', models[["tidy"]][["Exponential"]][["estimate"]][2], 'x', ')')
              ),
              #反函数
              ifelse(
                models[["tidy"]][["Inverse"]][["estimate"]][2] >= 0,
                paste0('y=', models[["tidy"]][["Inverse"]][["estimate"]][1],
                       '+', models[["tidy"]][["Inverse"]][["estimate"]][2], '/x'),
                paste0('y=', models[["tidy"]][["Inverse"]][["estimate"]][1] ,
                       models[["tidy"]][["Inverse"]][["estimate"]][2], '/x')
              ),
              #幂函数
              paste0('y=', models[["tidy"]][["Power"]][["estimate"]][1] , '(',
                     'x^', models[["tidy"]][["Power"]][["estimate"]][2], ')')
            ),
            stringsAsFactors = FALSE
          )
          mod <- mod %>% left_join(df, by = c('model' = 'mod')) %>% dplyr::select(
            model,formula,p.value,everything())%>%
            mutate_if(is.numeric,  ~sprintf('%.3f',.))
          message('>>> These (', paste(df$mod, collapse = '-'), ') models are fitted')
          mod
        }
      }else if(all(x!=0,y>0)){
        modlist <- list(Linear = lm(y~x),#线性
                        Quadratic = lm(y~x+I(x^2)),#二次
                        Compound = lm(log(y)~x),#复合
                        Growth = lm(log(y)~x),#生长
                        S = lm(log(y)~I(1/x)),#S曲线
                        Exponential = lm(log(y)~x),#指数
                        Inverse = lm(y~I(1/x))#反函数
        )
        if(aug){
          models <- modlist %>% map(broom::augment)
        }else{
          models[['tidy']] <- modlist %>% map(broom::tidy)
          models[['glance']] <- modlist %>% map(broom::glance)
          #对系数进行还原
          models$tidy$Compound[,2] <- exp(models$tidy$Compound[,2])#复合
          models$tidy$Exponential[1,2] <- exp(models$tidy$Exponential[1,2])#指数
          #生成表达式过程-----
          #1系数小数点保留,为提取方程的系数做准备
          models[["tidy"]] <- models[["tidy"]] %>% map(function(df){
            round(df %>% dplyr::select_if(is.numeric),2)
          })
          #提取模型为数据框
          models$glance %>% bind_rows(.id = "model") -> mod
          #生成方程
          df <- data.frame(
            mod =  c(
              "Linear",
              "Quadratic",
              "Compound",
              "Growth",
              "S",
              "Exponential",
              "Inverse"
            ),
            formula = c(
              #线性
              ifelse(
                models[["tidy"]][["Linear"]][["estimate"]][2] >= 0,
                paste0('y=', models[["tidy"]][["Linear"]][["estimate"]][1], '+',
                       models[["tidy"]][["Linear"]][["estimate"]][2], 'x'),
                paste0('y=', models[["tidy"]][["Linear"]][["estimate"]][1] ,
                       models[["tidy"]][["Linear"]][["estimate"]][2], 'x')
              ),
              #二次
              if (all(models[["tidy"]][["Quadratic"]][["estimate"]][2:3] >=
                      0)) {
                paste0('y=', models[["tidy"]][["Quadratic"]][["estimate"]][1],
                       '+', models[["tidy"]][["Quadratic"]][["estimate"]][2], 'x',
                       '+', models[["tidy"]][["Quadratic"]][["estimate"]][3], 'x^2')
              } else if (models[["tidy"]][["Quadratic"]][["estimate"]][2] >=
                         0 &
                         models[["tidy"]][["Quadratic"]][["estimate"]][3] <=
                         0) {
                paste0('y=', models[["tidy"]][["Quadratic"]][["estimate"]][1], '+',
                       models[["tidy"]][["Quadratic"]][["estimate"]][2], 'x',
                       models[["tidy"]][["Quadratic"]][["estimate"]][3], 'x^2')
              } else if (models[["tidy"]][["Quadratic"]][["estimate"]][2] <=
                         0 &
                         models[["tidy"]][["Quadratic"]][["estimate"]][3] >=
                         0) {
                paste0('y=', models[["tidy"]][["Quadratic"]][["estimate"]][1],
                       models[["tidy"]][["Quadratic"]][["estimate"]][2], 'x', '+',
                       models[["tidy"]][["Quadratic"]][["estimate"]][3], 'x^2')
              } else if (models[["tidy"]][["Quadratic"]][["estimate"]][2] <=
                         0 &
                         models[["tidy"]][["Quadratic"]][["estimate"]][3] <=
                         0) {
                paste0('y=', models[["tidy"]][["Quadratic"]][["estimate"]][1],
                       models[["tidy"]][["Quadratic"]][["estimate"]][2], 'x',
                       models[["tidy"]][["Quadratic"]][["estimate"]][3], 'x^2')
              },
              #复合
              paste0('y=', models[["tidy"]][["Compound"]][["estimate"]][1] , '(',
                     models[["tidy"]][["Compound"]][["estimate"]][2], '^x', ')'),
              #生长
              ifelse(
                models[["tidy"]][["Growth"]][["estimate"]][2] >= 0,
                paste0('y=', 'exp(', models[["tidy"]][["Growth"]][["estimate"]][1] ,
                       '+', models[["tidy"]][["Growth"]][["estimate"]][2], 'x', ')'),
                paste0('y=', 'exp(', models[["tidy"]][["Growth"]][["estimate"]][1] ,
                       models[["tidy"]][["Growth"]][["estimate"]][2], 'x', ')')
              ),
              #S
              ifelse(
                models[["tidy"]][["S"]][["estimate"]][2] >= 0,
                paste0('y=', 'exp(', models[["tidy"]][["S"]][["estimate"]][1] , '+',
                       models[["tidy"]][["S"]][["estimate"]][2], '/x', ')'),
                paste0('y=', 'exp(', models[["tidy"]][["S"]][["estimate"]][1] ,
                       models[["tidy"]][["S"]][["estimate"]][2], '/x', ')')
              ),
              #指数函数
              ifelse(
                models[["tidy"]][["Exponential"]][["estimate"]][2] >= 0,
                paste0('y=', models[["tidy"]][["Exponential"]][["estimate"]][1] ,
                       'exp(', models[["tidy"]][["Exponential"]][["estimate"]][2], 'x',
                       ')'),
                paste0('y=', models[["tidy"]][["Exponential"]][["estimate"]][1] ,
                       'exp(', models[["tidy"]][["Exponential"]][["estimate"]][2], 'x', ')')
              ),
              #反函数
              ifelse(
                models[["tidy"]][["Inverse"]][["estimate"]][2] >= 0,
                paste0('y=', models[["tidy"]][["Inverse"]][["estimate"]][1],
                       '+', models[["tidy"]][["Inverse"]][["estimate"]][2], '/x'),
                paste0('y=', models[["tidy"]][["Inverse"]][["estimate"]][1] ,
                       models[["tidy"]][["Inverse"]][["estimate"]][2], '/x')
              )),
            stringsAsFactors = FALSE
          )
          mod <- mod %>% left_join(df, by = c('model' = 'mod')) %>% dplyr::select(
            model,formula,p.value,everything())%>%
            mutate_if(is.numeric,  ~sprintf('%.3f',.))
          message('>>> These (', paste(df$mod, collapse = '-'), ') models are fitted')
          mod
        }
      }else if(any(x<=0, y<=0)){
        modlist <- list(Linear = lm(y~x),#线性
                        Quadratic = lm(y~x+I(x^2))#二次
        )
        if(aug){
          models <- modlist %>% map(broom::augment)
        }else{
          models[['tidy']] <- modlist %>% map(broom::tidy)
          models[['glance']] <- modlist %>% map(broom::glance)
          #生成表达式过程-----
          #1系数小数点保留,为提取方程的系数做准备
          models[["tidy"]] <- models[["tidy"]] %>% map(function(df){
            round(df %>% dplyr::select_if(is.numeric),2)
          })
          #提取模型为数据框
          models$glance %>% bind_rows(.id = "model") -> mod
          #生成方程
          df <- data.frame(
            mod =  c(
              "Linear",
              "Quadratic"
            ),
            formula = c(
              #线性
              ifelse(
                models[["tidy"]][["Linear"]][["estimate"]][2] >= 0,
                paste0('y=', models[["tidy"]][["Linear"]][["estimate"]][1], '+',
                       models[["tidy"]][["Linear"]][["estimate"]][2], 'x'),
                paste0('y=', models[["tidy"]][["Linear"]][["estimate"]][1] ,
                       models[["tidy"]][["Linear"]][["estimate"]][2], 'x')
              ),
              #二次
              if (all(models[["tidy"]][["Quadratic"]][["estimate"]][2:3] >=
                      0)) {
                paste0('y=', models[["tidy"]][["Quadratic"]][["estimate"]][1],
                       '+', models[["tidy"]][["Quadratic"]][["estimate"]][2], 'x',
                       '+', models[["tidy"]][["Quadratic"]][["estimate"]][3], 'x^2')
              } else if (models[["tidy"]][["Quadratic"]][["estimate"]][2] >=
                         0 &
                         models[["tidy"]][["Quadratic"]][["estimate"]][3] <=
                         0) {
                paste0('y=', models[["tidy"]][["Quadratic"]][["estimate"]][1], '+',
                       models[["tidy"]][["Quadratic"]][["estimate"]][2], 'x',
                       models[["tidy"]][["Quadratic"]][["estimate"]][3], 'x^2')
              } else if (models[["tidy"]][["Quadratic"]][["estimate"]][2] <=
                         0 &
                         models[["tidy"]][["Quadratic"]][["estimate"]][3] >=
                         0) {
                paste0('y=', models[["tidy"]][["Quadratic"]][["estimate"]][1],
                       models[["tidy"]][["Quadratic"]][["estimate"]][2], 'x', '+',
                       models[["tidy"]][["Quadratic"]][["estimate"]][3], 'x^2')
              } else if (models[["tidy"]][["Quadratic"]][["estimate"]][2] <=
                         0 &
                         models[["tidy"]][["Quadratic"]][["estimate"]][3] <=
                         0) {
                paste0('y=', models[["tidy"]][["Quadratic"]][["estimate"]][1],
                       models[["tidy"]][["Quadratic"]][["estimate"]][2], 'x',
                       models[["tidy"]][["Quadratic"]][["estimate"]][3], 'x^2')
              }
            ),
            stringsAsFactors = FALSE
          )
          mod <- mod %>% left_join(df, by = c('model' = 'mod')) %>% dplyr::select(
            model,formula,p.value,everything())%>%
            mutate_if(is.numeric,  ~sprintf('%.3f',.))
          message('>>> These (', paste(df$mod, collapse = '-'), ') models are fitted')
          mod
        }
      }else{
        message('\n>>>Please Check for both zero and negative values in the data')
      }
    }else{
      #进行logistic分布回归
      #自启动模型_提供参数----
      k <- nls(y ~ SSlogis(x, Asym, xmid, scal), data = df) %>%
        coef() %>% .[[1]]#k
      a <-
        1 / (nls(y ~ SSlogis(x, Asym, xmid, scal), data = df) %>%
               coef() %>% .[[3]])#a
      b <- nls(y ~ SSlogis(x, Asym, xmid, scal), data = df) %>%
        coef() %>% .[[2]]#b
      #生成方程
      models[['formula']] <- data.frame(
        mod = "Logistic",
        formula = #"Logistic分布"
          ifelse(
            b >= 0,
            paste0('y=', k, '/(', '1+exp(-', a, '(x-', b, ')))'),
            paste0('y=', k, '/(', '1+exp(-', a, '(x+', -b, ')))')
          ),
        stringsAsFactors = FALSE
      )
      #获得参数值
      log_mod <- nls(y ~ k / (1 + exp(-a * (x - b))),
                     start = list(k = k, a = a, b = b))
      models[['tidy']] <- log_mod %>% broom::tidy()
      models[['glance']] <- log_mod %>% broom::glance() %>%
        cbind(.,
              #计算r2，非线性回归的r2的适用性有待商榷
              data.frame(r.squared_f =
                           (1 - sum((
                             y - fitted(log_mod)
                           ) ^ 2) / sum((
                             y - mean(y)
                           ) ^ 2))))
      if (!is.null(models)) {
        message('The Logistic model is fitted')
      }
      models
    }
  }else{
    stop('>>>Please output numeric variables(both x and y)')
  }
}

#' Show the formula list
#' @description Show the formula list
#' @return show_formulalist(): Returns a data.frame
#' @examples show_formulalist()
#' @export
show_formulalist <- function() {
  mod_list = data.frame(
    mod_names = c(
      "Linear",
      "Quadratic",
      "Compound",
      "Growth",
      "Logarithmic",
      "S",
      "Exponential",
      "Inverse",
      "Power",
      "Logistic"
    ),
    CN_names = c(
      "线性",
      "二次",
      "复合",
      "生长",
      "对数",
      "S",
      "指数",
      "反函数",
      "幂函数",
      "Logistic分布"
    ),
    formula = c(
      'y=b0+b1x',
      'y=b0+b1x+b2x^2' ,
      'y=b0(b1^x)',
      'y=e^(b0+b1x)',
      'y=b0+b1ln(x)',
      'y=e^(b0+b1/x)',
      'y=b0e^(b1x)',
      'y=b0+b1/x',
      'y=b0(x^b1)',
      'y=k/(1+exp(-a(x-b)))'
    )
  )
  print(mod_list)
}
