gini <- function(x) {
    n <- length(x)
    xx <- sort(x)
    2 * (sum(seq(1, n) * xx)) / (n * sum(xx)) - 1
}

#function omnivory(A; weighted = true)
    ## Convert to Bool if preference matrix
    #tlvl = trophic_levels(A .!= 0)
    #omnivory = []
    #for i in 1:size(A, 1)
        #link_indexes = findall(!=(0), A[i,:])
        #prey_tlvl = tlvl[link_indexes]
        ## Relative preference:
        #if sum(A[i,:]) == 0 # if no prey
            #rel_pref = []
            #push!(omnivory, 0)
            #out = 0
        #else
            ## To vector if sparse vector
            #prey_tlvl_var = (prey_tlvl .- mean(prey_tlvl)).^2
            #if weighted
                #rel_pref = Vector(A[i,link_indexes] ./ sum(A[i,:]))
                #out = sum(prey_tlvl_var .* rel_pref)
            #else
                #out = sum(prey_tlvl_var)
            #end
            #push!(omnivory, out)
        #end
    #end
    #omnivory
#end
