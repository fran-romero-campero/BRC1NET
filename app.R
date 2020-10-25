## Author:  Francisco J. Romero-Campero
## Contact: Francisco J. Romero-Campero - fran@us.es 

## Input to test 
## input <- list(selected.multiple.tfs = c("BRC1", "HEC1", "ABI5"), cluster = "1")

## options(repos = BiocInstaller::biocinstallRepos())
## getOption("repos")

# Run the application 
shinyApp(ui = ui, server = server)

