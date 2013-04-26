#! /usr/bin/env ruby
#
# Filter for FASTA files
#

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

      opts.on("--codonize",
              "Trim sequence to be at multiple of 3 nucleotides") do |b|
        options.codonize = b
      end

      opts.on("--min [size]",
              "Set minimum sequence size") do |min|
        options.min = min.to_i
      end

      opts.on("-v", "--[no-]verbose", "Run verbosely") do |v|
        options.verbose = v
      end

      opts.separator ""
      opts.separator "Other options:"

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

FastaReader::emit_fastarecord(-> { ARGF.gets }) { | rec |
  if options.codonize
    size = rec.seq.size
    rec.seq = rec.seq[0..size - (size % 3) - 1]
  end
  next if options.min and req.seq.size < options.min
  
  print rec.to_fasta
}

