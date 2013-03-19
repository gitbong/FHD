FUNCTION healpix_combine_inds,hpx_cnv,hpx_inds=hpx_inds,reverse_ind_reference=reverse_ind_reference
;assumes hpx_cnv is a pointer array to hpx_cnv structures generated by healpix_cnv_generate.pro

n_obs=N_Elements(hpx_cnv)
hpx_min=Lon64arr(n_obs) 
hpx_max=Lon64arr(n_obs) 
FOR obs_i=0,n_obs-1 DO BEGIN
    hpx_min[obs_i]=Min((*hpx_cnv[obs_i]).inds)
    hpx_max[obs_i]=Max((*hpx_cnv[obs_i]).inds)
ENDFOR
hpx_min=Min(hpx_min)
hpx_max=Max(hpx_max)
n_hpx=hpx_max-hpx_min+1
n_hpx_full=nside2npix((*hpx_cnv[obs_i]).nside)

ind_hist=lonarr(n_hpx)
hist_arr=Ptrarr(n_obs)
FOR obs_i=0,n_obs-1 DO BEGIN
    ind_hist1=histogram((*hpx_cnv[obs_i]).inds,min=hpx_min,max=hpx_max,/bin)
    hist_arr[obs_i]=Ptr_new(ind_hist1)
    ind_hist+=ind_hist1
ENDFOR
ind_hist1=0 ;free memory

ind_use=where(ind_hist,n_hpx_use)
hpx_inds=ind_use+hpx_min

hpx_ind_map=Ptrarr(size(hpx_cnv,/dimension))
FOR obs_i=0,n_obs-1 DO BEGIN
    ind_use2=where((*hist_arr[obs_i])[ind_use],n_use2)
    hpx_ind_map[obs_i]=Ptr_new(ind_use2)
    *hist_arr[obs_i]=0  ;free memory
ENDFOR
Ptr_free,hist_arr

IF Arg_Present(reverse_ind_reference) THEN BEGIN
    reverse_ind_reference=Lon64arr(n_hpx_full)-1
    reverse_ind_reference[hpx_inds]=L64indgen(N_Elements(hpx_inds))
;    full_ind_reference=Lon64arr(hpx_max-hpx_min+1)-1
;    full_ind_reference[ind_use]=ind_use
ENDIF
RETURN,hpx_ind_map
END