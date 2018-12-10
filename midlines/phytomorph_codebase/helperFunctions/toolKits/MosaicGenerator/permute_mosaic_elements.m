function mosaic_el_perm=permute_mosaic_elements(needed_no_of_elements,no_of_existing_elms)

mosaic_el_perm=1:no_of_existing_elms;

if needed_no_of_elements>no_of_existing_elms
    mosaic_el_perm=imresize(mosaic_el_perm,[1,needed_no_of_elements],'nearest');
end
mosaic_el_perm=mosaic_el_perm(randperm(length(mosaic_el_perm)));
mosaic_el_perm=mosaic_el_perm(1:needed_no_of_elements);