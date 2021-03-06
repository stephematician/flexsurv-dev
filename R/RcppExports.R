# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

dgenf_work <- function(x, mu, sigma, Q, P, log) {
    .Call('flexsurv_dgenf_work', PACKAGE = 'flexsurv', x, mu, sigma, Q, P, log)
}

pgenf_work <- function(q, mu, sigma, Q, P, lower_tail, give_log) {
    .Call('flexsurv_pgenf_work', PACKAGE = 'flexsurv', q, mu, sigma, Q, P, lower_tail, give_log)
}

check.genf <- function(mu, sigma, Q, P) {
    .Call('flexsurv_check_genf', PACKAGE = 'flexsurv', mu, sigma, Q, P)
}

dgengamma_work <- function(x, mu, sigma, Q, log) {
    .Call('flexsurv_dgengamma_work', PACKAGE = 'flexsurv', x, mu, sigma, Q, log)
}

pgengamma_work <- function(q, mu, sigma, Q, lower_tail, give_log) {
    .Call('flexsurv_pgengamma_work', PACKAGE = 'flexsurv', q, mu, sigma, Q, lower_tail, give_log)
}

check.gengamma <- function(mu, sigma, Q) {
    .Call('flexsurv_check_gengamma', PACKAGE = 'flexsurv', mu, sigma, Q)
}

dgompertz_work <- function(x, shape, rate, log) {
    .Call('flexsurv_dgompertz_work', PACKAGE = 'flexsurv', x, shape, rate, log)
}

pgompertz_work <- function(q, shape, rate, lower_tail, give_log) {
    .Call('flexsurv_pgompertz_work', PACKAGE = 'flexsurv', q, shape, rate, lower_tail, give_log)
}

check.gompertz <- function(shape, rate) {
    .Call('flexsurv_check_gompertz', PACKAGE = 'flexsurv', shape, rate)
}

dllogis_work <- function(x, shape, scale, log) {
    .Call('flexsurv_dllogis_work', PACKAGE = 'flexsurv', x, shape, scale, log)
}

pllogis_work <- function(q, shape, scale, lower_tail, give_log) {
    .Call('flexsurv_pllogis_work', PACKAGE = 'flexsurv', q, shape, scale, lower_tail, give_log)
}

check.llogis <- function(shape, scale) {
    .Call('flexsurv_check_llogis', PACKAGE = 'flexsurv', shape, scale)
}

exph <- function(y) {
    .Call('flexsurv_exph', PACKAGE = 'flexsurv', y)
}

dexph <- function(y) {
    .Call('flexsurv_dexph', PACKAGE = 'flexsurv', y)
}

basis_vector <- function(knots, x) {
    .Call('flexsurv_basis_vector', PACKAGE = 'flexsurv', knots, x)
}

basis_matrix <- function(knots, x) {
    .Call('flexsurv_basis_matrix', PACKAGE = 'flexsurv', knots, x)
}

dbasis_vector <- function(knots, x) {
    .Call('flexsurv_dbasis_vector', PACKAGE = 'flexsurv', knots, x)
}

dbasis_matrix <- function(knots, x) {
    .Call('flexsurv_dbasis_matrix', PACKAGE = 'flexsurv', knots, x)
}

