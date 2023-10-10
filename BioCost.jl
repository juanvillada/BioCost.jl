using BioSequences
using CSV
using DataFrames
using FASTX
using Statistics
using Transducers
using BangBang

@time begin

println("")
println("📭\tStarting BioCost.jl")
println("🚀\tUsing $(Threads.nthreads()) thread(s)")

include("data/cost_dictionaries/dictionaryaminoacids.jl")
include("data/cost_dictionaries/dictionarydna.jl")

path_to_dir_with_genomes = ARGS[1]
output_file = ARGS[2]
sequence_type = ARGS[3]

println("📂\tAnalyzing folder: $path_to_dir_with_genomes")

vector_of_genome_paths = readdir( abspath(path_to_dir_with_genomes), join=true)
println("🗃\tThe folder contains $(length(vector_of_genome_paths)) genomes")


function compute_AA_synthesis_cost(my_fasta_seq::FASTX.FASTA.Record)
    cost_C = 0
    cost_H = 0
    cost_O = 0
    cost_N = 0
    cost_S = 0

    seq_id = identifier(my_fasta_seq)
    # println(seq_id)
    my_seq = sequence(LongAA, my_fasta_seq) # The first arg makes sure it opens as AA Type

    for each_aa in my_seq
        cost_C += get(dict_AA_C , each_aa, 0)
        cost_H += get(dict_AA_H , each_aa, 0)
        cost_O += get(dict_AA_O , each_aa, 0)
        cost_N += get(dict_AA_N , each_aa, 0)
        cost_S += get(dict_AA_S , each_aa, 0)
    end

    length_unAmb = n_certain(my_seq)
    length_Amb   = n_ambiguous(my_seq)
    length_total = length_unAmb + length_Amb
    cost_C = cost_C / length_unAmb
    cost_H = cost_H / length_unAmb
    cost_O = cost_O / length_unAmb
    cost_N = cost_N / length_unAmb
    cost_S = cost_S / length_unAmb

    return Dict(
        "seq_id" => seq_id,
        "length_total" => length_total,
        "length_unAmb" => length_unAmb,
        "length_Amb"   => length_Amb,
        "C_per_AA" => cost_C,
        "H_per_AA" => cost_H,
        "O_per_AA" => cost_O,
        "N_per_AA" => cost_N,
        "S_per_AA" => cost_S
    )
end


function compute_nucleotide_synthesis_cost(my_fasta_seq::FASTX.FASTA.Record)
    cost_C = 0
    cost_H = 0
    cost_O = 0
    cost_N = 0
    cost_S = 0

    seq_id = identifier(my_fasta_seq)
    my_seq = sequence(LongDNA{4}, my_fasta_seq)

    for each_nt in my_seq
        cost_C += get(dict_DNA_C , each_nt, 0)
        cost_H += get(dict_DNA_H , each_nt, 0)
        cost_O += get(dict_DNA_O , each_nt, 0)
        cost_N += get(dict_DNA_N , each_nt, 0)
    end

    length_unAmb = n_certain(my_seq)
    length_Amb   = n_ambiguous(my_seq)
    length_total = length_unAmb + length_Amb
    cost_C = cost_C / length_unAmb
    cost_H = cost_H / length_unAmb
    cost_O = cost_O / length_unAmb
    cost_N = cost_N / length_unAmb

    return Dict(
        "seq_id" => seq_id,
        "length_total" => length_total,
        "length_unAmb" => length_unAmb,
        "length_Amb"   => length_Amb,
        "C_per_nt" => cost_C,
        "H_per_nt" => cost_H,
        "O_per_nt" => cost_O,
        "N_per_nt" => cost_N
    )
end

function compute_cost_for_a_genome(path_to_genome::String)
    genome_name = splitext(basename(path_to_genome))[1]
    # println(genome_name)

    cost_df = DataFrame()

    n_seqs = 0

    # open(FASTA.Reader, path_to_genome) do fasta_file
    #     for each_seq in fasta_file
    #         append!(
    #             cost_df,
    #             DataFrame( compute_AA_synthesis_cost(each_seq) )
    #         )
    #         n_seqs += 1
    #     end
    # end

    open(FASTA.Reader, path_to_genome) do fasta_file
        each_seq = FASTA.Record()
        while !eof(fasta_file)
            read!(fasta_file, each_seq)
            if sequence_type == "AA"
                append!(
                    cost_df,
                    DataFrame( compute_AA_synthesis_cost(each_seq) )
                )
            elseif sequence_type == "DNA"
                append!(
                    cost_df,
                    DataFrame( compute_nucleotide_synthesis_cost(each_seq) )
                )
            end
            n_seqs += 1
        end
    end

    if sequence_type == "AA"
        columns_to_process = [:length_total, :length_Amb, :length_unAmb, :C_per_AA, :H_per_AA, :O_per_AA, :N_per_AA, :S_per_AA]
    elseif sequence_type == "DNA"
        columns_to_process = [:length_total, :length_Amb, :length_unAmb, :C_per_nt, :H_per_nt, :O_per_nt, :N_per_nt]
    end

    cost_df = combine(
        cost_df, columns_to_process .=> mean,
        columns_to_process .=> std
    )

    # cost_df = combine(
    #     cost_df, columns_to_process .=> mean(filter(!isnan, x)),
    #     columns_to_process .=> std(filter(!isnan, x))
    # )

    insertcols!(cost_df, 1, :genome_name => genome_name)
    insertcols!(cost_df, 2, :n_seqs => n_seqs)

    return cost_df
end

final_cost_df = foldxt(append!!, Map(compute_cost_for_a_genome), vector_of_genome_paths)

# final_cost_df = DataFrame()
# for each_genome in vector_of_genome_paths
#     append!(final_cost_df, compute_cost_for_a_genome( each_genome ))
# end

CSV.write(output_file, final_cost_df)

println("📬\tOutput file saved to: $(output_file)")

println("✅ \tDone!")

println("")

end # timing end

println("")


# TODO:
# multiple dispatch for FILE vs DIR
# other amino acids (e.g. X)
# is X counted by n_ambiguous? Yes!
# What about gaps? "-" is not counted either by n_certain() nor
# n_ambiguous() in AA seqs
