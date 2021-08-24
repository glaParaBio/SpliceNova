gff_attribute <- function(attribute_str, key) {
    xkey <- paste0('^', key, '=')
    lst <- strsplit(attribute_str, ';')

    value <- rep('', length(lst))
    for(i in 1:length(lst)) {
        kv <- lst[[i]][grepl(xkey, lst[[i]])]
        stopifnot(length(kv) <= 1)
        if(length(kv) == 1) {
            value[i] <- sub(xkey, '', kv)
        } 
    }
    return(value)
}
