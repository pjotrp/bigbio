#! /usr/bin/env ruby
#
# Filter for FASTA files
#

require 'bigbio'

count = 0
FastaReader::emit_fastarecord(-> { gets }) { | rec |
  p count += 1
  print rec.to_fasta
}

