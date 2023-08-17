getGeneID <- function(RegionList, Path_To_AnnoSqlite, Path_To_OrgdbSqlite){
  ### RegionList: a list of regions
  ### Path_To_AnnoSqlite: path to obtain annotation file
  ### Path_To_OrgdbSqlite: path to obtain mappings between different version of gene names
  ### this function is used to annotate which gene the peak is located on
  
  txdb = loadDb(Path_To_AnnoSqlite)
  Orgdb = loadDb(Path_To_OrgdbSqlite)

  allgene = genes(txdb, columns="gene_id", single.strand.genes.only=TRUE)
  features <- AnnotationDbi::select(Orgdb, keys = allgene$gene_id, 
                                    columns = c("SYMBOL"), 
                                    keytype = "ENTREZID") 

  genes.GR = GRanges(seqnames = Rle(allgene@seqnames), allgene@ranges, strand = Rle(allgene@strand),
                     gene_id = Rle(allgene$gene_id), geneSymbol = Rle(features$SYMBOL))
  
  regions.GR = GRanges(seqnames = Rle(RegionList$chr), IRanges(RegionList$start, RegionList$end),
                       strand = Rle(RegionList$strand))
  hits = findOverlaps(regions.GR, genes.GR)
  query = queryHits(hits)
  subject = subjectHits(hits)
  ENTREZID = NULL
  Symbol = NULL
  for(i in 1:length(regions.GR)){
    idx = which(query==i)
    if (length(idx) == 0){
      ENTREZID[i] = NA
      Symbol[i] = NA
    }else{
      ENTREZID[i] = paste0(genes.GR$gene_id[subject[idx]], collapse = " / ")
      Symbol[i] = paste0(genes.GR$geneSymbol[subject[idx]], collapse = " / ")
    }
  }
  
  res = data.frame(chr = regions.GR@seqnames, start = regions.GR@ranges@start, 
                   end = regions.GR@ranges@start + regions.GR@ranges@width -1,
                   strand = regions.GR@strand,
                   gene_Symbol = Symbol
                   #gene_ENTREZID = ENTREZID
                   )
  return(res)
}

