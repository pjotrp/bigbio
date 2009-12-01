# Ruby DocTest
#
# Test getorf functionality.
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

BigBio's ORF predictor uses the rapid Biolib::Emboss::transeq function to
translate Nucleotide sequences to Amino Acid sequences. Next it allows
several heuristics to select potential ORF's and returns them with their
reading frame and position in the nucleotide sequence. 

One of the advantages of PredictORF is that it returns both the amino acid and
matching nucleotide sequences.

  >> require 'bigbio/sequence/predictorf'

  >> id = "PUT-157a-Arabidopsis_thaliana-126"
  >> descr = "PlantGDB Arabidopsis_thaliana Jan_15_2007"
  >> sequence = "ATCATTAGCAACACCAGCTTCCTCTCTCTCGCTTCAAAGTTCACTACTCGTGGATCTCGT
                 CTTCAGTGTACAGTATCAAGGGCTCGATCTGCGGTGGATGAGACATCAGATTCAGGAGCT
                 TTTCAAAGAACTGCATCGACATCCGTAACTTCGTTTCAAAAGATTCCAATTCTCAGTTTC
                 AGCTGAATCTGGTAGATACCATCTTTACATATCGTATGCTTGTCATGGGCTTCTAGATGC
                 CTTTCATACTTAAAGATCAAAGGACTTGACGATGCAATAAGCTTCTCGTCTGTAAAACCC"

Pick up the EMBOSS translation table:

  >> trn_table = Biolib::Emboss.ajTrnNewI(1)

Initiate the ORF prediction class

  >> predict = PredictORF.new(id,descr,sequence,trn_table)
 
The methods return a list of ORFCandidate ordered by sequence length. Here
we look for all ORF's between STOP codons (with a minimal size of 30 AA).

  >> orflist = predict.stopstop
  >> orflist.size
  => 9

Get the first (and largest) ORF

  >> orf = orflist.first

The id contains the number of the ORF at the last position (like EMBOSS' 
getorf does)

  >> orf.id
  => "PUT-157a-Arabidopsis_thaliana-126_1"

The description contains 'XX' for the STOPSTOP search. Unlike getorf it shows
the reading frame.

  >> orf.descr
  => "[XX 0 - 183; +1] PlantGDB Arabidopsis_thaliana Jan_15_2007"
  >> orf.aa.seq.size
  => 61
  >> orf.aa.seq
  => "IISNTSFLSLASKFTTRGSRLQCTVSRARSAVDETSDSGAFQRTASTSVTSFQKIPILSFS"

The ORF object contains more information:

  >> orf.nt.start
  => 0
  >> orf.frame
  => 1

The matching sequence with the AA

  >> orf.nt.seq
  => "ATCATTAGCAACACCAGCTTCCTCTCTCTCGCTTCAAAGTTCACTACTCGTGGATCTCGTCTTCAGTGTACAGTATCAAGGGCTCGATCTGCGGTGGATGAGACATCAGATTCAGGAGCTTTTCAAAGAACTGCATCGACATCCGTAACTTCGTTTCAAAAGATTCCAATTCTCAGTTTCAGCT"

And it keeps track of the full nucleotide sequence

  >> orf.nt.fullseq
  => "ATCATTAGCAACACCAGCTTCCTCTCTCTCGCTTCAAAGTTCACTACTCGTGGATCTCGTCTTCAGTGTACAGTATCAAGGGCTCGATCTGCGGTGGATGAGACATCAGATTCAGGAGCTTTTCAAAGAACTGCATCGACATCCGTAACTTCGTTTCAAAAGATTCCAATTCTCAGTTTCAGCTGAATCTGGTAGATACCATCTTTACATATCGTATGCTTGTCATGGGCTTCTAGATGCCTTTCATACTTAAAGATCAAAGGACTTGACGATGCAATAAGCTTCTCGTCTGTAAAACCC"

Let's check one of the others

  >> orf = orflist[3]
  >> orf.nt.start
  => 99
  >> orf.nt.stop
  => 231
  >> orf.frame
  => 3
  >> orf.aa.seq
  => "DIRFRSFSKNCIDIRNFVSKDSNSQFQLNLVDTIFTYRMLVMGF"
  >> orf.nt.seq
  => "GAGACATCAGATTCAGGAGCTTTTCAAAGAACTGCATCGACATCCGTAACTTCGTTTCAAAAGATTCCAATTCTCAGTTTCAGCTGAATCTGGTAGATACCATCTTTACATATCGTATGCTTGTCATGGGCTT"
  >> orf.nt.fullseq[orf.nt.start..orf.nt.stop]
  => "GAGACATCAGATTCAGGAGCTTTTCAAAGAACTGCATCGACATCCGTAACTTCGTTTCAAAAGATTCCAATTCTCAGTTTCAGCTGAATCTGGTAGATACCATCTTTACATATCGTATGCTTGTCATGGGCTT"

Naming for each ORF

  >> orf.id
  => "PUT-157a-Arabidopsis_thaliana-126_3"
  >> orf.descr
  => "[XX 99 - 231; +3] PlantGDB Arabidopsis_thaliana Jan_15_2007"

If you want to know the size of the smallest ORF

  >> orflist.last.aa.seq.size
  => 30

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
