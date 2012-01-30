# Ruby DocTest
#
# Translate a nucleotide sequence into six reading frames using the fast
# EMBOSS transeq function
#
# Run with ./runner.rb or
#
#   env DATA=../../../biolib/src/test/data/emboss/ ruby ../../../biolib/tools/rubydoctest/bin/rubydoctest test_getorf.rb
#
# or
#
#   ruby ../../../biolib/tools/rubydoctest/bin/rubydoctest test_getorf.rb
#
# Documentation with rd2 -r rd/rd2html-lib *.rb

cwd = File.dirname(__FILE__)
Dir.chdir(cwd)

$: << "../../lib"

if $UNITTEST

=begin

BibBio's frame translation uses the rapid Biolib::Emboss::transeq function to
translate Nucleotide sequences to Amino Acid sequences. 

  >> require 'bigbio/sequence/translate'

  >> id = "PUT-157a-Arabidopsis_thaliana-126"
  >> descr = "PlantGDB Arabidopsis_thaliana Jan_15_2007"
  >> sequence = "ATCATTAGCAACACCAGCTTCCTCTCTCTCGCTTCAAAGTTCACTACTCGTGGATCTCGT
                 CTTCAGTGTACAGTATCAAGGGCTCGATCTGCGGTGGATGAGACATCAGATTCAGGAGCT
                 TTTCAAAGAACTGCATCGACATCCGTAACTTCGTTTCAAAAGATTCCAATTCTCAGTTTC
                 AGCTGAATCTGGTAGATACCATCTTTACATATCGTATGCTTGTCATGGGCTTCTAGATGC
                 CTTTCATACTTAAAGATCAAAGGACTTGACGATGCAATAAGCTTCTCGTCTGTAAAACCC"


   >> translate = Nucleotide::Translate.new(1)
   >> list = translate.aa_frames(sequence)

We should have six frames

   >> list.size
   => 6

   >> aa = list.first
   >> aa[:frame]
   => 1

This result matches the one from the EMBOSS web interface:

   >> aa[:sequence] 
   => "IISNTSFLSLASKFTTRGSRLQCTVSRARSAVDETSDSGAFQRTASTSVTSFQKIPILSFS*IW*IPSLHIVCLSWASRCLSYLKIKGLDDAISFSSVKP"

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
