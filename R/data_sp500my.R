#' Daily log return of the GSPC from 2004-01-05 to 2024-08-02
#'
#' GSPC is a ticker for S&P 500 Index, see \url{https://en.wikipedia.org/wiki/S\%26P_500#cite_note-10}.
#'
#' @name SP500MY
#' @docType data
#' @usage data(SP500MY)
#'
#' @format A xts object with length 5180 of daily data for 1 variable.
#' \describe{
#'   \item{log return}{log return}
#'   ...
#' }
#'
#' @keywords datasets
#' @source \url{https://finance.yahoo.com/quote/\%5EGSPC/history?period1=949363200&period2=1025222400&interval=1d&filter=history&frequency=1d}
#' @importFrom xts xts
"SP500MY"
