# Ruby DocTest
#
# Test paired NA+AA sequence functionality.
#
# Run with ./runner.rb or
#
#   ../../../../tools/rubydoctest/bin/rubydoctest test_paired.rb
#
# Documentation with rd2 -r rd/rd2html-lib *.rb

cwd = File.dirname(__FILE__)
Dir.chdir(cwd)

# $: << '../../../mappings/swig/ruby/rqtl/'

require 'biolib/biolib_core'
Biolib::Biolib_core.biolib_setloglevel(7)

if $UNITTEST

=begin

  >> $: << '..'
  >> require 'bio/sequence2'

Sequence pairs are paired NA and AA sequences where one can be tested against
the other (through translation) and an NA sequence can be aligned against an AA
alignment (protein alignment to nucleotide alignment, also known as pal2nal).

=end


$: << '..'
require 'db/fasta'
require 'test/unit'

TESTDIR = '../../../test/data/fasta'
NA_FILE = TESTDIR + "/na.fa"
AA_FILE = TESTDIR + "/aa.fa"

class TestBiolibFasta < Test::Unit::TestCase

  def setup
  end

  def test_indexer
    na_in = FastaReader.new(NA_FILE, :regex => '(\d+)\s', :index => true)
    rec = na_in.get("122")
    assert_equal("122",rec.id)
    assert_equal("121",na_in.get("121").id)
  end

end

end # $UNITTEST
