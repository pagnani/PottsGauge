using Documenter, PottsGauge

push!(LOAD_PATH,"../src/")
makedocs(sitename = "PottsGauge",
	 modules = [PottsGauge],
	 doctest = true)
deploydocs(
	   branch = "gh-pages",
	   repo = "github.com/pagnani/PottsGauge.git",
	   versions = ["stable" => "v^"]
	  )
