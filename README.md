# BIGBIO

BigBio = BIG DATA for Ruby

BigBio is an initiative to a create high performance low-memory
libraries for big data computing in biology.

BigBio may use BioLib C/C++/D functions for increasing performance and
reducing memory consumption.

In a way, this is an experimental project. I use it for
experimentation, but what is in here should work fine. If you wish to
contribute subscribe to the BioRuby and/or BioLib mailing lists
instead.

# Overview

* BigBio can translate nucleotide sequences to amino acid
  sequences using an EMBOSS C function, or BioRuby's translator.
* BigBio has a terrific FASTA file emitter which iterates FASTA files and
  iterates sequences without loading everything in memory. There is
  also an indexed edition
* BioBio has a flexible FASTA filter 
* BigBio has an ORF emitter which parses DNA/RNA sequences and emits
  ORFs between START_STOP or STOP_STOP codons.
* BigBio has a Phylip (PAML style) emitter and writer

# Installation

The easy way

```sh
gem install bio-bigbio
```

in your code

```ruby
require 'bigbio'
```

# Command line tools

Some functionality comes also as executable command line tools (see the
./bin directory). Use the -h switch to get information. Current tools
are 

1. getorf: fetch all areas between start-stop and stop-stop codons in six frames (using EMBOSS when biolib is available)
2. nt2aa.rb: translate in six frames (using EMBOSS when biolib is available)
3. fasta_filter.rb

## Command line Fasta Filter

The CLI filter accepts standard Ruby commands. 

Filter sequences that contain more than 25% C's

```sh
fasta_filter.rb --filter "rec.seq.count('C') > rec.seq.size*0.25" test/data/fasta/nt.fa
```

Look for IDs containing -126 and sequences ending on CCC

```sh
fasta_filter.rb --filter "rec.id =~ /-126/ or rec.seq =~ /CCC$/" test/data/fasta/nt.fa
```

Filter out all masked sequences that contain more than 10% masked
nucleotides

```sh
fasta_filter.rb --filter "rec.seq.count('N')<rec.seq.size*0.10" 
```

Next to rec.id and rec.seq, you have rec.descr and 'num' as variables,
so to skip every other record

```sh
fasta_filter.rb --filter "num % 2 == 0" 
```

Find all sequences that contain a stop codon in the sequence

```sh
fasta_filter.rb --filter 'rec.seq =~ /\*./' aa.fa
```

Rewrite all sequences to lower case, you can use the useful rewrite
option

```sh
fasta_filter.rb --rewrite 'rec.seq = rec.seq.downcase'
```

Rewrite the FASTA descriptors

```sh
fasta_filter.rb --rewrite 'rec.descr =~ /gene=(\S+)/; rec.descr = $1' test.fa
```

Filters and rewrites can be combined. The rest is up to your
imagination! One final example to remove low quality sequences from an
amino acid file (one amino acid dominates):

```sh
fasta_filter.rb --filter "count = {} ; rec.seq.each_char { |c| count[c] ||= 1 ; count[c] += 1 }; count.values.max.to_f/rec.seq.size < 0.30" < aa.fa > aa1.fa
```


# API Examples

## Iterate through a FASTA file

Read a file without loading the whole thing in memory

```ruby
require 'bigbio'

fasta = FastaReader.new(fn)
fasta.each do | rec |
  print rec.descr,rec.seq
end
```

Since FastaReader parses the ID, write a tab file with id and sequence

```ruby
i = 1
print "num\tid\tseq\n"
FastaReader.new(fn).each do | rec |
  if rec.id =~ /(AT\w+)/
    print i,"\t",$1,"\t",rec.seq,"\n"
    i += 1
  end
end
```

wich, for example, can be turned into RDF with the
[bio-table](https://github.com/pjotrp/bioruby-table) biogem.

## Write a FASTA file

Write a FASTA file. The simple way

```ruby
fasta = FastaWriter.new(fn)
fasta.write('Test',"agtcta")
```

Any object can be passed in, however, as long as it responds to
'descr' and 'seq.to_s', or 'id' and 'seq.to_s'. E.g.

```ruby
class StorageObject
  attr_accessor :descr, :seq
end

mysequence = StorageObject.new
mysequence.descr = 'Test'
mysequence.seq = "agtcta"
```

write the FASTA file

```ruby
fasta = FastaWriter.new(fn)
fasta.write(mysequence)
```

## Transform a FASTA file

You can combine above FastaReader and FastaWriter to transform
sequences, e.g.

```ruby
fasta = FastaWriter.new(in_fn)
FastaReader.new(out_fn).each do | rec |
  # Strip the description down to the second ID
  (id1,id2) = /(\S+)\s+(\S+)/.match(rec.descr)
  fasta.write(id2,rec.seq)
end
```

The downside to this approach is the explicit file naming. What if you
want to use STDIN or some other source instead? I have come round to
the idea of using a combination of lambda and block. For example:

```ruby
  FastaReader::emit_fastarecord(-> {gets}) { |rec|
    print FastaWriter.to_fasta(rec)
  }
```

which takes STDIN line by line, and outputs FASTA on STDOUT. This is 
a better design as the FastaReader and FastaWriter know nothing of
the mechanism fetching and displaying data. These can both be 'pure' 
functions. Note also that the data is never fully loaded into RAM.

Here the transformer functional style

```ruby
  FastaReader::emit_fastarecord(-> {gets}) { |rec|
    (id1,id2) = /(\S+)\s+(\S+)/.match(rec.descr)
    print FastaWriter.to_fasta(id2,req.seq)
  }
```

## Fetch ORFs from a sequence

BigBio can parse a sequence for ORFs. Together with the FastaReader
little memory gets used

```ruby
predictorf = PredictORF.new(id,descr,"ATCATTAGCAACACCAGCTTCCTCTCTCTCGCTTCAAAGTTCACTACTCGTGGATCTCGT")
# get all ORFs between start and stop codons, longer than 30 bps
orfs = predictorf.startstop(30)
# get all sequences between stop codons
seqs = predictorf.stopstop(0)
```

## Rapid DNA/RNA to amino acid translation

Translate with EMBOSS C library, if linked, otherwise use BioRuby

```ruby
trn_table = Bio::Big::TranslationAdapter.translation_table(1)
translate = Nucleotide::Translate.new(trn_table)
aa_frames = translate.aa_6_frames("ATCATTAGCAACACCAGCTTCCTCTCTCTCGCTTCAAAGTTCACTACTCGTGGATCTCGT")
```

## Walk a FASTA (reference) genome

Genomes and BACS often come as large (continuous) FASTA files. When
variant/position queries happen on sorted data, the genome can be
walked through once reading the whole file serially. This is what
FastaGenomeReader does.

The following code assumes the FASTA descriptors contain 

```ruby
  >13 dna:chromosome chromosome:GRCh37:13:1:115169878:1
```

so 'chr' is captured, as well as 'start' and 'stop'. Using 
[bio-vcf](https://github.com/pjotrp/bioruby-vcf):

```ruby
genome = FastaGenomeReader.new('Hs_GRCh37_gatk.fasta', -> 
  { |descr| a = skip,skip,skip,chr,start,stop = descr.split(':')
      chr, start.to_i, stop.to_i } )

STDIN.each_line do | line |
  next if line =~ /^#/
  fields = VcfLine.parse(line)
  rec = VcfRecord.new(fields,header)
  if rec.var == genome.ref(vcf.chr,vcf.pos+1)
    # do something
  end
end
```

FastaGenomeReader is buffered and tiled. You can override the size of
64K.

# Project home page

Information on the source tree, documentation, examples, issues and
how to contribute, see

  http://github.com/pjotrp/bigbio

The BioRuby community is on IRC server: irc.freenode.org, channel: #bioruby.

# Cite

If you use this software, please cite one of
  
* [BioRuby: bioinformatics software for the Ruby programming language](http://dx.doi.org/10.1093/bioinformatics/btq475)
* [Biogem: an effective tool-based approach for scaling up open source software development in bioinformatics](http://dx.doi.org/10.1093/bioinformatics/bts080)

# Biogems.info

This Biogem is published on [biogems.info](http://biogems.info/index.html#bigbio)

# Copyright

Copyright (c) 2011-2014 Pjotr Prins. See LICENSE for further details.

