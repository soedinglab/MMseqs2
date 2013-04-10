#!/usr/bin/ruby

if ARGV[0] == "-h" || ARGV.size != 3
    puts "\nGenerates fake meaningsless prefiltering scores lists in the right format.\n\nUSAGE: create_pref_scores_dummy.rb uniprot_formatted_input_file.fas num_pref_seqs_per_seq outdir\n"
    exit
end

infile = ARGV[0].dup
n = ARGV[1].dup.to_i
outdir = ARGV[2].dup

maxscore = 100

puts "Reading database..."
inf = File.new(infile, 'r')
# read the UniProt IDs
ids = Array.new
while (line=inf.gets)
    next if line[0,1] != '>'
    ids << line[4,6]
end
inf.close

if n > ids.size
    puts "Too small db for the given prefiltering list length."
    exit
end

puts "Generating random prefiltering lists..."
# members of the prefiltering lists, first empty
members = Hash.new
ids.each { |id|
    members[id] = Array.new
}

ids.each_with_index { |id, i|
    cand = ids[i+1,ids.length-i]
    cand.shuffle!

    id_n = members[id].size
    cand.each { |m|
        next if members[m].size == n
        break if members[id].size == n
        score = (rand*maxscore*100).round/100.0
        members[id] << [m, score]
        members[m] << [id, score]
    }
}

puts "Sorting prefiltering lists..."
# sort prefiltering lists by the score
members.each_key { |k|
    members[k].sort! { |a,b| b[1] <=> a[1] }
}


members.keys.each { |id|
    next if members[id].size == 0
    outpath = outdir + "/" + id
    outf = File.new(outpath, 'w')
    id_mems = members[id]
    id_mems.each{ |m|
        outf << m[0]
        outf << " "
        outf << m[1]
        outf << "\n"
    }
    outf.close
}

