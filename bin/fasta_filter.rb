#! /usr/bin/env ruby
#
# Filter for FASTA files
#

$: << File.dirname(__FILE__)+'/../lib'

require 'bigbio'
require 'optparse'
require 'ostruct'

class OptParser

  #
  # Return a structure describing the options.
  #
  def self.parse(args)
    # The options specified on the command line will be collected in *options*.
    # We set default values here.
    options = OpenStruct.new
    options.codonize = false
    options.verbose = false

    opt_parser = OptionParser.new do |opts|
      opts.banner = "Usage: fasta_filter.rb [options]"

      opts.separator ""
      opts.separator "Specific options:"

      opts.on("--filter expression","Filter on Ruby expression") do |expr|
        options.filter = expr
      end

      opts.on("--codonize",
              "Trim sequence to be at multiple of 3 nucleotides") do |b|
        options.codonize = b
      end

      opts.on("--min size",
              "Set minimum sequence size") do |min|
        options.min = min.to_i
      end

      opts.on("--id","Write out ID only") do |b|
        options.id = b
      end

      opts.on("-v", "--[no-]verbose", "Run verbosely") do |v|
        options.verbose = v
      end

      opts.separator ""
      opts.separator "Examples:"
      opts.separator ""
      opts.separator "  fasta_filter.rb --filter \"rec.id =~ /-126/ or rec.seq =~ /CCC$/\" test/data/fasta/nt.fa"
      opts.separator "  fasta_filter.rb --filter \"rec.seq.count('C') > rec.seq.size*0.25\" test/data/fasta/nt.fa"
      opts.separator "  fasta_filter.rb --filter \"rec.descr =~ /C. elegans/\" test/data/fasta/nt.fa"
      opts.separator "  fasta_filter.rb --filter \"num % 2 == 0\" test/data/fasta/nt.fa"
      opts.separator ""
      opts.separator "Other options:"
      opts.separator ""

      opts.on_tail("-h", "--help", "Show this message") do
        puts opts
        exit
      end

    end

    opt_parser.parse!(args)
    options
  end  # parse()

end  # class OptParser

options = OptParser.parse(ARGV)

num = -1
FastaReader::emit_fastarecord(-> { ARGF.gets }) { | rec |
  num += 1
  if options.filter
    b = eval(options.filter)
    next if not eval(options.filter)
  end
  if options.codonize
    size = rec.seq.size
    rec.seq = rec.seq[0..size - (size % 3) - 1]
  end
  next if options.min and rec.seq.size < options.min
  rec.descr = rec.id if options.id
  
  print rec.to_fasta
}

