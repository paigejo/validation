# script for counting downloads of CRAN packages
# month = last 30 days
getDownloads = function(packageNames, type=c('total', 'month', 'week', 'day')) {
  type = match.arg(type)
  
  if(FALSE) {
    require(devtools)
    devtools::install_github("metacran/cranlogs")
  }
  require(cranlogs)
  
  if(type == "total") {
    sum(cran_downloads(packages = packageNames, from = "1990-06-30")$count)
  } else if(type == "month") {
    sum(cran_downloads(when = "last-month", packages = packageNames)$count)
  } else if(type == "week") {
    sum(cran_downloads(when = "last-week", packages = packageNames)$count)
  } else if(type == "day") {
    sum(cran_downloads(when = "last-day", packages = packageNames)$count)
  }
}