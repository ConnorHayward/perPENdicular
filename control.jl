abs_range = 5:2:100
alpha_range = 0:0.02:1
n_events = 100000

for abs in abs_range
    for α in alpha_range
    script = "/PEN/det/setDetName perPENdicularr_test
       /PEN/det/setABS $abs
       /PEN/det/setSigAlpha $α
       /run/beamOn $n_events"

       open("test.mac","w") do f
           write(f,script)
       end
        run(`./PEN -m test.mac`)
    end
end

using Glob, DataFrames, Plots, CSV
folder = "/home/iwsatlas1/hayward/Documents/perPENdicular_sim/Data/"
files = glob("*.csv",folder);

all_abs = Float64[]
all_alphas = Float64[]
all_ratio = Float64[]

matrix_pairs = Dict{Array{Float64},Float64}()

for file in values(files)
    ratio = Float64
#    try
        data = CSV.read(file,datarow=7)

#        println(length(data[1]))
       ratio = length(data[:1]) / n_events
   # catch UnDefVarError
   #     ratio = 0

    file = replace(file,","=>".")
    start_ind = findnext("_",file[length(folder):end],1)[1]+length(folder);

    tmp_string = file[start_ind:end];

    start_ind = start_ind+findnext("_",tmp_string[1:end],1)[1];
    abs = parse(Float64,file[start_ind:start_ind+5])
    start_ind = start_ind+findnext("_",tmp_string[1:end],1)[1];
    start_ind = start_ind+findnext("_",tmp_string[1:end],1)[1];
    start_ind = start_ind+findnext("_",tmp_string[1:end],1)[1];

    if abs == 100
        alpha = parse(Float64,file[start_ind-1:start_ind+3])
    elseif abs == 5
        alpha = parse(Float64,file[start_ind-3:start_ind+3])
    else
        alpha = parse(Float64,file[start_ind-2:start_ind+3])
    end

    push!(all_abs,abs)
    push!(all_alphas,alpha)
    push!(all_ratio,ratio)

    push!(matrix_pairs,[abs, alpha]=>ratio)
end

filter_abs = sort(unique(all_abs))
filter_alphas = sort(unique(all_alphas))

matrix = Matrix(undef,length(filter_abs),length(filter_alphas))

for value_abs in values(filter_abs)
    for value_alphas in values(filter_alphas)
        row  = findfirst(x->x==value_abs, filter_abs)
        col  = findfirst(x->x==value_alphas, filter_alphas)

        matrix[row,col]=matrix_pairs[[filter_abs[row],filter_alphas[col]]]

    end
end
savefig(heatmap(filter_abs,filter_alphas,matrix,ylabel="Polish Paramater, alpha[A.U]",xlabel="ABS Length [mm]",ylims=[0,1],xlims=[5:100]),"Documents/perPENdicular_sim/matrix_plot.pdf")
