#!/usr/bin/env ruby
#
# fasta_sort: Sorts a FASTA file and outputs sorted unique records as FASTA again
#
# Usage:
#
#   fasta_sort inputfile(s)

require 'bio'

include Bio

table = Hash.new  
ARGV.each do | fn |
  Bio::FlatFile.auto(fn).each do | seq |
    table[seq.definition] ||= seq.data
  end
end

table.sort.each do | definition, data |
  rec = Bio::FastaFormat.new('> '+definition.strip+"\n"+data)
  print rec
end

