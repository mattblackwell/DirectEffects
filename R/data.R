#' Data on historical plough use and the socioeconomic status of
#' women.
#'
#' A dataset to replicate the analysis in Alesina, Giuliano, and Nunn
#' (2013).
#'
#' \itemize{
#'   \item isocode. 3-letter code for the country.
#' \item flfp2000. Female labor force participation in 2000
#' \item female_ownership. Percent of firms with female ownership (in
#'   latest survey year)
#' \item women_politics. Women in Politics in 2000, WDI
#' \item plow. Animal plow cultivation variable (v39): Using
#'   Ethnologue - pop weighted
#' \item agricultural_suitability. overall (millets, sorghum, wheat,
#'   barley, rye): share defined as suitable
#' \item tropical_climate. Frac land: tropics and subtropics: using
#'   Ethnologue - pop weighted
#' \item large_animals. presence of large animals
#' \item political_hierarchies. Jurisdictional hierarchy beyond local
#'   community (v33): Using Ethnologue - pop weighted
#' \item economic_complexity. Settlement patterns (v30)
#' \item ln_income. ln (income)
#' \item ln_income_squared. ln (income) ^2
#' \item country. country name
#' \item communist_dummy. Communism indicator variable
#' \item rugged. Ruggedness (Terrain Ruggedness Index, 100 m.)
#' \item years_interstate_conflict. Years of interstate conflict, 1800-2007 - from COW
#' \item serv_va_gdp2000. Value Added in Service/GDP in 2000
#' \item polity2_2000. Polity 2 measure taken from the Polity IV dataset
#' \item oil_pc. oil production/GDP
#' \item ... other variables as annotated in the source.
#' }
#'
#' @docType data
#' @name ploughs
#' @usage data(ploughs)
#' @format A data frame with 234 observations and 57 variables.
#' @source \url{http://qje.oxfordjournals.org/content/128/2/469}
#' @references Alesina, A., Giuliano, P., & Nunn, N. (2013). On the
#' Origins of Gender Roles: Women and the Plough. The Quarterly
#' Journal of Economics, 128(2), 469-530.
NULL