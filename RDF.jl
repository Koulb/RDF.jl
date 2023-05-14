using Plots

include("routines.jl")

"""
This script calculates and plots g(r) using the analytic
normalization procedure
"""

#DASHBOARD

RMax = 100 # maximal radial distance
DeltaR = 0.2 # bin width
NBins = 1000 # number of equally spaced bins
InOutDir = "test/"
InitailData = InOutDir * "blobs_expBSA_19_Position_1_SIRT20i_1n.inv.csv"
ParticleFile = InOutDir * "test.dat"
SimuName = ParticleFile[1:end-4]

function readPartCoor(PartFile)
	# open the .txt file and read in all datalines
	DataFile = open(PartFile, "r")
	datalines = readlines(DataFile)
	close(DataFile)

	# exctract all particle coordinates
    len = size(datalines)[1]
	Particles = Matrix{Float64}(undef, len, 3)

	for (index, line) in enumerate(datalines)
	   temp = split(line, " ")
	   Part = [parse(Float64,temp[1]) parse(Float64,temp[2]) parse(Float64,temp[3])]
       Particles[index,: ] = Part

	end
    
	return Particles
end

function plotRDF(r, Gr, SimuName)
    plot(r, Gr, color="green", label="", xlabel="r", ylabel="g(r)")
    plot!(r, ones(length(r)), color="black", linestyle=:dash, linewidth=1, label="")
    savefig(SimuName * "_Gr.pdf")
    return true
end

function saveRDFdata(rList, Gr, DeltaR, SimuName, FilePath)
	# create a new file for the raw data
	DataFile = open(SimuName * "_RDF.txt", "w")

	# write the header
	write(DataFile, "Radial Distribution Function for Experiment: ", SimuName, "\n")
	write(DataFile, "Original Path: ", FilePath, "\n")
	write(DataFile, "Radial Bin Width: ", string(DeltaR), "\n")
	write(DataFile, "Radial Distribution Function:\n")
	write(DataFile, "#r [nm]\t g(r)\n")

	# write the data
	for k in eachindex(rList)

	    write(DataFile, string(rList[k]), "\t", string(Gr[k]), "\n")

	end

	close(DataFile)

	return true
end

function ProcessPartCoor(PartFile, rList, DeltaR, SimuName)
	println("Current FilePath: ", PartFile)

	Particles = readPartCoor(PartFile)

	# calculate the RDF
	Gr = RDF_AnalyticNorm(Particles, rList, DeltaR)

	# save the RDF data to a file
	saveRDFdata(rList, Gr, DeltaR, SimuName, PartFile)

	# plot RDF
	plotRDF(rList, Gr, SimuName)

	return true
end

# Perform rotation of data 
println("Performing rotation of data...")
args = InitailData * " " * ParticleFile
run(`python rotate_axes.py $args`)
println("Rotation complete")

#create the list with radial distances
rList = range(0, stop=RMax, length=NBins)

println("Computing RDF:")
ProcessPartCoor(ParticleFile, rList, DeltaR, SimuName)

println("Program Finished")





