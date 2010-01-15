# Indexer module for the FASTA class
#
# This is a simple memory based key storage
#

module Indexer

  # Start using the indexer
  def indexer_use state
    if state
      @indexer = {}
    end
  end

  def indexer_set key, fpos
    raise "Trying to use 'set' when there is no index" if @indexer == nil
    raise "Indexer key #{key} alread in use for <#{@indexer[key]}>!" if @indexer[key]
    # p [key, fpos]
    @indexer[key] = fpos
  end

  # Get the key, return nil when not found
  def indexer_get key
    raise "Trying to use 'get' when there is no index" if @indexer == nil
    # raise "Indexer key #{key} not found!" if !@indexer[key]
    @indexer[key] 
  end

  def indexer_get_by_index idx
    @indexer.sort {|a,b| a[1]<=>b[1]} [idx]
  end
end

