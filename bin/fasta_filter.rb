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

      opts.on("--rewrite expression","Rewrite expression") do |expr|
        options.rewrite = expr
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
      opts.separator "  fasta_filter.rb test/data/fasta/nt.fa --rewrite 'rec.seq.downcase!'"
      opts.separator "  fasta_filter.rb --rewrite 'rec.descr =~ /gene=(\S+)/; rec.descr = $1' test.fa"
      opts.separator "  fasta_filter.rb --filter 'rec.seq =~ /\*./' aa.fa"
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
  # --- Filtering
  next if options.filter and not eval(options.filter)
  if options.codonize 
    # --- Round sequence to nearest 3 nucleotides
    size = rec.seq.size
    rec.seq = rec.seq[0..size - (size % 3) - 1]
  end
  # --- Only use sequences from MIN size
  next if options.min and rec.seq.size < options.min
  # --- Truncate description to ID
  rec.descr = rec.id if options.id
 
  # --- rewrite
  eval(options.rewrite) if options.rewrite
  print rec.to_fasta
}

