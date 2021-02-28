# reads DataFrame from baseFile and appendFile
# appends 2nd file to the 1st and writes result to baseFile
# So (warning) you might want to make a copy of the base file
# MGP Feb 2021


function PS_appendDataFrame(baseFile::String, appendFile::String)

    # read DataFrames
    DFbase = CSV.read(baseFile, DataFrame)
    N = DFbase.rep[end]
    println("Read ", N, " rows from ", baseFile)

    DFappend = CSV.read(appendFile, DataFrame)
    N2 = DFappend.rep[end]
    println("Read ", N2, " rows from ", appendFile)

    # increment rep column for combined reps
    DFappend.rep[:] .+= N

    # append
    DFnew = reduce(vcat, [DFbase, DFappend])

    # write combined data to base file
    CSV.write(baseFile, DFnew)
    println("Wrote ", N+N2, " rows to ", baseFile)

end

