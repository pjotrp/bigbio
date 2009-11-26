#! /usr/bin/ruby
#
# Performance testing routing for translating a FASTA file into six
# reading frames using the Bioruby routines.
#
# by Pjotr Prins (c) 2009
# 
# To reduce the impact of file IO you can run multiple iterations using
# a command line switch. 
#
# Usage: 
#
#   ruby -Ipath_to_bioruby translate_with_bioruby.rb --iter 10 nucleotides.fasta
#
# Example:
#
#   time ruby -I~/izip/git/opensource/bioruby/lib/ translate_with_bioruby.rb --iter 10 ../data/fasta/nt.fa > test.out
#
# Renders on my machine:
#
#  96 records 5760 times translated!
#   real    0m6.414s
#   user    0m5.928s
#   sys     0m0.384s


$: << '../../lib'

require 'bio'
require 'bigbio'

iter=1
fn = ARGV.shift

if fn == '--iter'
  iter = ARGV.shift.to_i
  fn = ARGV.shift
end

nt = FastaReader.new(fn)

nt.each { | rec |
  (0..iter).each do | repeat |
    seq = Bio::Sequence::NA.new(rec.seq)
    [-3,-2,-1,1,2,3].each do | frame |
      print "> ",rec.id," ",frame.to_s,"\n"
      print seq.translate(frame),"\n"
    end
  end
}

$stderr.print nt.size," records ",nt.size*6*iter," times translated!"




