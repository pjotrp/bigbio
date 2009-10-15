# Ruby DocTest
#
# Test paired nt+AA sequence functiontlity.
#
# Run with ./runner.rb or
#
#   ../../../../tools/rubydoctest/bin/rubydoctest test_paired.rb
#
# Documentation with rd2 -r rd/rd2html-lib *.rb

cwd = File.dirname(__FILE__)
Dir.chdir(cwd)

# $: << '../../../mappings/swig/ruby/rqtl/'

# require 'biolib/biolib_core'
# Biolib::Biolib_core.biolib_setloglevel(7)

if $UNITTEST

=begin

  >> $: << '..'
  !> require 'bio/sequence2'

Sequence pairs are paired NT and AA sequences where one can be tested against
the other (through translation) and an nt sequence can be aligned against an AA
alignment (protein alignment to nucleotide alignment, also known as pal2ntl).

=end


$: << '..'
require 'db/fasta'
require 'test/unit'

TESTDIR = '../../../test/data/fasta'
nt_FILE = TESTDIR + "/nt.fa"
AA_FILE = TESTDIR + "/aa.fa"

class TestBiolibFasta < Test::Unit::TestCase

  def setup
  end

  def test_indexer
    nt_in = FastaReader.new(nt_FILE, :regex => '(\d+)\s', :index => true)
    rec = nt_in.get("122")
    assert_equal("122",rec.id)
    assert_equal("121",nt_in.get("121").id)
  end

end

end # $UNITTEST
