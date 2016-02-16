user_guide.md : user_guide.Rmd
	# rm -rf cache/
	R --slave -e 'knitr::knit("user_guide.Rmd")'
	# R --slave -e 'rmarkdown::render("user_guide.Rmd")'
	perl -i -pe 's{figures/}{https://github.com/raychaudhurilab/custard/blob/master/vignettes/figures/}' user_guide.md
