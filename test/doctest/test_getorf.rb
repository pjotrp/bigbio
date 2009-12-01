# Ruby DocTest
#
# Test getorf functionality.
#
# Run with ./runner.rb or
#
#   ruby ../../../biolib/tools/rubydoctest/bin/rubydoctest test_getorf.rb
#
# Documentation with rd2 -r rd/rd2html-lib *.rb

cwd = File.dirname(__FILE__)
Dir.chdir(cwd)

if $UNITTEST

=begin

BigBio's ORF predictor uses the rapid Biolib::Emboss::transeq function to
translate Nucleotide sequences to Amino Acid sequences. Next it allows
several heuristics to select potential ORF's and returns them with their
reading frame and position in the nucleotide sequence. 

One of the advantages of PredictORF is that it returns both the amino acid and
nucleotide sequences.

  >> require 'bigbio'
  >> predict = PredictORF.new(sequence,table)
 
The methods return a list of ORFCandidate ordered by sequence length. Here
we look for all ORF's between STOP codons (with a minimal size of 30 AA).

  >> orflist = predict.stopstop
  >> orflist.size
  => 12

Get the first (and largest) ORF

  >> orf = orflist.first
  >> orf.name
  => "xxxxx"
  >> orf.id
  => "122"
  >> orf.aa.size
  => 100
  >> orf.aa.seq
  => "AAAA"

The ORF object contains more information:

  >> orf.nt.start
  => 34
  >> orf.frame
  => -1
  >> orf.nt.seq
  => "actg"
  >> orf.nt.fullseq
  => "actxxxx"

If you want to know the size of the smallest ORF

  >> orflist.last.aa.size
  => 33

=end


$: << '..'
require 'db/fasta'
require 'test/unit'

TESTDIR = '../../../test/data/fasta'
AA_FILE = TESTDIR + "/aa.fa"

# class TestBiolibFasta < Test::Unit::TestCase
# 
#   def setup
#   end
# 
#   def test_indexer
#   end
# 
# end

end # $UNITTEST
