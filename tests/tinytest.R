
if ( requireNamespace("tinytest", quietly=TRUE) ){
  home <- length(unclass(packageVersion("TruncatedNormal"))[[1]]) == 4
  tinytest::test_package("TruncatedNormal", at_home = home)
}

