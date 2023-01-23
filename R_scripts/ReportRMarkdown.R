# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve; Tiffeny De Coninck
# Maintainer: <limet@ugent.be>
# Script: ReportRMarkdown


##########Configuration##########

#Working directory

## @knitr NAAM
## options
PATH_USER <- 'xxx/Pipeline_metabolomics' # example path user

CODE_AUTORUN <- 'run code in ternminal automatically'
CODE_DEVELOPMENT <- 'run code manually in Rstudio for development locally'

#Adjustments
PATH <- PATH_USER
CODE_RUN_MODE <- CODE_AUTORUN


if (CODE_RUN_MODE == CODE_AUTORUN){
  #Recognize projectname from Rscript commando
  name_project = commandArgs(trailingOnly=TRUE) #name is given in "Rscript main.r projectname" as argument
  # test if there is one argument: if not, return an error
  if (length(args) == 0) {
    stop("Your projectname must be supplied", call.=FALSE)
  } 
  if (length(args) > 1) {
    stop("Only one projectname must be supplied", call.=FALSE)
  }
}
if (CODE_RUN_MODE == CODE_DEVELOPMENT){
  #give projectname in configuration.r (for code in development @Rstudio)
  name_project <- EXPERIMENT #'20200326_LA_REIMS_Avecom_pos9'
}

#Source input and output 
path_data_in <- file.path(PATH, 'Data/Input', name_project) #directory must exist!
path_data_out <- file.path(PATH, 'Data/Output', name_project) #moet reeds aangemaakt zijn (dus pipeline gelopen)


## @knitr STOP
##########Report RMarkdown##########
print(Sys.time())
start_time_total <- Sys.time()
print("Report RMarkdown - start!")


#Reportgeneration do NOT adjust

path_R_scripts <- file.path(PATH, 'R_scripts')
#path_R_scripts <- '/home/lca/lib/gedit'
dir.create(file.path(PATH, 'Data/Reports', name_project))

path_data_report <- file.path(PATH, 'Data/Reports', name_project)

Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio/bin/pandoc")

Configuration <- paste(PATH, '/Data/Input/' , name_project, '/Configuration.R', sep="")
source(Configuration)

if(RUN_PART_UNIVARIATE == RUN_CODE){
  rmarkdown::render(file.path(path_R_scripts, 'Report_univariate_analysis.Rmd'), output_file = file.path(path_data_report, 'Report_univariate_analysis'), output_format="html_document")
}

print("Report RMarkdown - done!")
print(Sys.time())
end_time_total <- Sys.time()
print(end_time_total - start_time_total)

#
####################
