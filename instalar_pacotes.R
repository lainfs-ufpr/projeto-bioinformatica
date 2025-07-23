instalar_pacotes <- function(pkgs, installer) {
  
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      installer(pkg)
    }
  }
}
