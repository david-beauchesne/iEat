tanimoto <- function(taxon1, taxon2) {
    # Tanimoto similarity measure, which compares two vectors x and y with n = |x| = |y| elements, and is defined as the size of the intersection (∩) of two sets divided by their union (∪):
    # The order of vectors taxon1 or taxon2 has no importance, as long as elements in vectors are unique
    if(length(taxon1) == 0 || length(taxon2) == 0) {   # If either length of taxon1 or taxon2 == 0, similarity == 0
        return(0.0)
    } else if(taxon1 == "" || taxon2 == "") { # "" or NAs, to verify
        return(0.0)
    } else {
        inter <- length(taxon2[match(taxon1, taxon2, nomatch = 0)])
        return(inter / ((length(taxon1) + length(taxon2)) - inter))
    }#if
}#end tanimoto function
