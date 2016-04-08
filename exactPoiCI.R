exactPoiCI <- function (x, conf.level=0.95) 
{
  poisson.test(x, conf.level = conf.level)$conf.int[1:2]
}