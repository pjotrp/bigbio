#! /usr/bin/env ruby
#
# Performance testing routing for translating a FASTA file into six
# reading frames using the Biolib (EMBOSS) routines.
#
# by Pjotr Prins (c) 2009
# 
# To reduce the impact of file IO you can run multiple iterations using
# a command line switch. 
#
# Usage: 
#
#   ruby -Ipath_to_biolib translate_with_biolib.rb --iter 10 nucleotides.fasta
#
# Example:
#
#   time ruby -I~/izip/git/opensource/biolib/lib/ translate_with_bioruby.rb --iter 10 ../data/fasta/nt.fa > test.out
#
# Renders on my machine:
#
#  96 records 5760 times translated!
#   real    0m0.290s
#   user    0m0.252s
#   sys     0m0.024s
#
# with a large file
#
#  22929 records 137574 times translated!
#   real    0m20.306s
#   user    0m15.997s
#   sys     0m1.344s


$: << '../../lib'

require 'biolib/emboss'
require 'bigbio'

iter=1
fn = ARGV.shift

if fn == '--iter'
  iter = ARGV.shift.to_i
  fn = ARGV.shift
end

nt = FastaReader.new(fn)
trnTable = Biolib::Emboss.ajTrnNewI(1);

nt.each { | rec |
  (0..iter).each do | repeat |
    ajpseq   = Biolib::Emboss.ajSeqNewNameC(rec.seq,"Test sequence")

    [-3,-2,-1,1,2,3].each do | frame |
      ajpseqt  = Biolib::Emboss.ajTrnSeqOrig(trnTable,ajpseq,frame)
      aa       = Biolib::Emboss.ajSeqGetSeqCopyC(ajpseqt)
      print "> ",rec.id," ",frame.to_s,"\n"
      print aa,"\n"
    end
  end
}

$stderr.print nt.size," records ",nt.size*6*iter," times translated!"




