def gff_attribute(attribute_str, key):
    # ID=exon_PBANKA_1307600.1-E1;Parent=PBANKA_1307600.1;gene_id=PBANKA_1307600;biotype=exon
    xkey = key + '='
    lst = attribute_str.split(';')
    kv = [x for x in lst if x.startswith(xkey)] 
    assert len(kv) <= 1
    if len(kv) == 0:
        return ''
    kv = kv[0]
    value = kv[len(xkey):]
    return value
