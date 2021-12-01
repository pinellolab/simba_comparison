
args = commandArgs(trailingOnly = TRUE)

input_h5ad = args[1]
metadata_file = args[2]
output_file = args[3]

if (length(args) < 3) stop("Usage: Rscript convert.R input_h5ad metadata_file output_file")

system(paste0("python convert_h5_matched.py ", input_h5ad, " .tmp_embedding.txt"))

dir.create(dirname(output_file))

barcode_idx = 3

tbl = read.table(".tmp_embedding.txt", head = T, sep = "\t", row.names =1)
if(metadata_file == "same"){ 
    convert_barcode = F
    barcode_idx = 1 } else {
    meta = read.table(metadata_file, head = T); 
    convert_barcode = T
}

barcodes = sapply(strsplit(rownames(tbl), split = "_|\\."), function(s) {paste0(s[barcode_idx:(length(s)-1)], collapse = ".")})
print(head(barcodes))
modalities = sapply(strsplit(rownames(tbl), split = "_"), function(s) s[2])


if (convert_barcode){
    stopifnot(barcodes %in% c(levels(meta$atac.bc), levels(meta$rna.bc)))

    rna_cell_id = sapply(barcodes, function(bc){ 
        if (bc %in% levels(meta$rna.bc)) {return(bc)} 
        else {
            rna.bc = levels(meta$rna.bc)[match(bc, levels(meta$atac.bc))]; 
            return(rna.bc) } })
} else{
    rna_cell_id = barcodes
}

stopifnot(table(rna_cell_id) == 2)

tbl$batch = rna_cell_id

modality_cell_id = sapply(1:length(barcodes), function(i){ 
                          bc = barcodes[i]
                          modality = modalities[i]
                          return(paste0(modality, ".", bc))
        }
)

save.image()
rownames(tbl) = modality_cell_id

write.table(tbl, file = output_file, sep = "\t", col.names = T, row.names = T, quote =F)
#system("rm .tmp_embedding.txt")
