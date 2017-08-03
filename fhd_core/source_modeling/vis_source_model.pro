FUNCTION vis_source_model,skymodel, obs, status_str, psf, params, vis_weight_ptr, cal, jones, model_uv_arr=model_uv_arr,$
    file_path_fhd=file_path_fhd, timing=timing, silent=silent, uv_mask=uv_mask, error=error, beam_arr=beam_arr,$
    fill_model_visibilities=fill_model_visibilities, use_pointing_center=use_pointing_center, vis_model_ptr=vis_model_ptr,$
    spectral_model_uv_arr=spectral_model_uv_arr,model_delay_filter=model_delay_filter, _Extra=extra
fill_model_visibilities=1
t0=Systime(1)
IF N_Elements(error) EQ 0 THEN error=0
IF N_Elements(file_path_fhd) EQ 0 THEN file_path_fhd=''
IF N_Elements(silent) EQ 0 THEN silent=1

IF N_Elements(skymodel) EQ 0 THEN fhd_save_io,status_str,skymodel,var='skymodel',/restore,file_path_fhd=file_path_fhd,_Extra=extra
IF N_Elements(obs) EQ 0 THEN fhd_save_io,status_str,obs,var='obs',/restore,file_path_fhd=file_path_fhd,_Extra=extra
IF N_Elements(psf) EQ 0 THEN fhd_save_io,status_str,psf,var='psf',/restore,file_path_fhd=file_path_fhd,_Extra=extra
IF N_Elements(params) EQ 0 THEN fhd_save_io,status_str,params,var='params',/restore,file_path_fhd=file_path_fhd,_Extra=extra
IF Min(Ptr_valid(vis_weight_ptr)) EQ 0 THEN fhd_save_io,status_str,vis_weight_ptr,var='vis_weights',/restore,file_path_fhd=file_path_fhd,_Extra=extra
IF N_Elements(jones) EQ 0 THEN fhd_save_io,status_str,jones,var='jones',/restore,file_path_fhd=file_path_fhd,_Extra=extra

IF Keyword_Set(skymodel) THEN BEGIN
    galaxy_flag=skymodel.galaxy_model
    diffuse_filepath=skymodel.diffuse_model
    n_sources=skymodel.n_sources
ENDIF ELSE BEGIN
    galaxy_flag=0
    diffuse_filepath=''
    n_sources=0
ENDELSE
heap_gc

pol_names=obs.pol_names

;extract information from the structures
n_pol=obs.n_pol
n_spectral=obs.degrid_spectral_terms
dimension=obs.dimension
elements=obs.elements
degpix=obs.degpix
kbinsize=obs.kpix
kx_span=kbinsize*dimension ;Units are # of wavelengths
ky_span=kx_span
icomp=Complex(0,1)

xvals=meshgrid(dimension,elements,1)-dimension/2
yvals=meshgrid(dimension,elements,2)-elements/2

IF Keyword_Set(uv_mask) THEN uv_mask_use=uv_mask ELSE uv_mask_use=Fltarr(dimension,elements)+1
uv_mask_use[*,elements/2+psf.dim:*]=0. 

freq_bin_i=(*obs.baseline_info).fbin_i

frequency_array=(*obs.baseline_info).freq

;Save the original structures, and setup variables for a double bandwidth model
if keyword_set(model_delay_filter) then begin
  obs_original = replicate(obs,1)
  params_original = replicate(params,1)
  freq_bin_i=[INTARR(obs.n_freq/2)+freq_bin_i[0],freq_bin_i,INTARR(obs.n_freq/2)+freq_bin_i[obs.n_freq-1]]
  for i=1,obs.n_freq/2 do frequency_array = [(*obs.baseline_info).freq[0]-obs.freq_res*i,frequency_array]
  for i=1,obs.n_freq/2 do frequency_array = [frequency_array,(*obs.baseline_info).freq[N_elements((*obs.baseline_info).freq)-1]+obs.freq_res*i]
endif

nfreq_bin=Max(freq_bin_i)+1
nbaselines=obs.nbaselines
n_samples=obs.n_time
n_freq=obs.n_freq

;Set up the structures such that a double bandwidth model is made
if keyword_set(model_delay_filter) then begin
  
  n_freq = n_freq*2
  freq_use = INTARR(n_freq) + 1
  
  ;Setup obs structure (switch to using the init pro?)
  arr={tile_A:(*obs.baseline_info).tile_a,tile_B:(*obs.baseline_info).tile_b,bin_offset:(*obs.baseline_info).bin_offset,Jdate:(*obs.baseline_info).Jdate,freq:Float(frequency_array),fbin_i:Long(freq_bin_i),$
    freq_use:Fix(freq_use),tile_use:(*obs.baseline_info).tile_use,time_use:(*obs.baseline_info).time_use,tile_names:(*obs.baseline_info).tile_names,tile_height:(*obs.baseline_info).tile_height,tile_flag:(*obs.baseline_info).tile_flag}
  obs={code_version:obs.code_version,instrument:obs.instrument,obsname:obs.obsname,$
    dimension:obs.dimension,elements:obs.elements,nbaselines:obs.nbaselines,dft_threshold:obs.dft_threshold,double_precision:obs.double_precision,$
    kpix:obs.kpix,degpix:obs.degpix,obsaz:obs.obsaz,obsalt:obs.obsalt,obsra:obs.obsra,obsdec:obs.obsdec,$
    zenra:obs.zenra,zendec:obs.zendec,obsx:obs.obsx,obsy:obs.obsy,zenx:obs.zenx,zeny:obs.zeny,$
    phasera:obs.phasera,phasedec:obs.phasedec,orig_phasera:obs.orig_phasera,orig_phasedec:obs.orig_phasedec,$
    n_pol:obs.n_pol,n_tile:obs.n_tile,n_tile_flag:obs.n_tile_flag,n_freq:Long(n_freq),n_freq_flag:0L,n_time:obs.n_time,n_time_flag:0L,$
    n_vis:obs.n_vis,n_vis_in:obs.n_vis_in,n_vis_raw:obs.n_vis_raw,nf_vis:obs.nf_vis,beam_integral:obs.beam_integral,beam2_integral:obs.beam2_integral,$
    pol_names:obs.pol_names,$
    jd0:obs.jd0,max_baseline:obs.max_baseline,min_baseline:obs.min_baseline,delays:obs.delays,lon:obs.lon,lat:obs.lat,alt:obs.alt,$
    freq_center:obs.freq_center,freq_res:obs.freq_res,time_res:obs.time_res,astr:obs.astr,alpha:obs.alpha,pflag:obs.pflag,cal:obs.cal,$
    residual:0,vis_noise:obs.vis_noise,baseline_info:Ptr_new(arr),meta_data:obs.meta_data,meta_hdr:obs.meta_hdr,$
    degrid_spectral_terms:obs.degrid_spectral_terms,grid_spectral_terms:obs.grid_spectral_terms,grid_info:obs.grid_info,healpix:obs.healpix}    
  
  psf=fhd_struct_init_psf(beam_ptr=psf.beam_ptr,fbin_i=freq_bin_i,xvals=psf.xvals,yvals=psf.yvals,$
    n_pol=psf.n_pol,n_freq=nfreq_bin,freq_cen=obs.freq_center,group_arr=psf.id, psf_dim=psf.dim, $
    psf_resolution = psf.resolution, complex_flag=psf.complex_flag)

endif

n_freq_bin=N_Elements(freq_bin_i)
IF N_Elements(vis_model_ptr) LT n_pol THEN vis_model_ptr=intarr(n_pol)
IF n_spectral EQ 0 THEN spectral_model_uv_arr=intarr(n_pol)

vis_dimension=Float(nbaselines*n_samples)

IF Min(Ptr_valid(model_uv_arr)) EQ 0 THEN BEGIN
    model_uv_arr=Ptrarr(n_pol,/allocate)
    FOR pol_i=0,n_pol-1 DO *model_uv_arr[pol_i]=Complexarr(dimension,elements)
ENDIF
IF (Min(Ptr_valid(spectral_model_uv_arr)) EQ 0) AND (n_spectral GT 0) THEN BEGIN
    spectral_model_uv_arr=Ptrarr(n_pol,n_spectral,/allocate)
    FOR pol_i=0,n_pol-1 DO FOR s_i=0,n_spectral-1 DO *spectral_model_uv_arr[pol_i,s_i]=Complexarr(dimension,elements)
ENDIF

IF n_sources GT 0 THEN BEGIN ;test that there are actual sources in the source list
    ;convert Stokes entries to instrumental polarization (weighted by one factor of the beam) 
    ;NOTE this is for record-keeping purposes, since the Stokes flux values will actually be used
    source_list=skymodel.source_list
    source_list.extend=Pointer_copy(source_list.extend)
    source_list=stokes_cnv(source_list,jones,beam_arr=beam_arr,/inverse,_Extra=extra) 
    model_uv_arr1=source_dft_model(obs,jones,source_list,t_model=t_model,sigma_threshold=2.,$
        spectral_model_uv_arr=spectral_model_uv_arr1,uv_mask=uv_mask_use,_Extra=extra)
    FOR pol_i=0,n_pol-1 DO *model_uv_arr[pol_i]+=*model_uv_arr1[pol_i];*uv_mask_use 
    FOR pol_i=0,n_pol-1 DO FOR s_i=0,n_spectral-1 DO *spectral_model_uv_arr[pol_i,s_i]+=*spectral_model_uv_arr1[pol_i,s_i];*uv_mask_use
    undefine_fhd,model_uv_arr1,spectral_model_uv_arr1,source_list
    IF ~Keyword_Set(silent) THEN print,"DFT timing: "+strn(t_model)+" (",strn(n_sources)+" sources)"
ENDIF


IF galaxy_flag THEN gal_model_uv=fhd_galaxy_model(obs,jones,spectral_model_uv_arr=gal_spectral_model_uv,antialias=1,/uv_return,_Extra=extra)
IF Min(Ptr_valid(gal_model_uv)) GT 0 THEN FOR pol_i=0,n_pol-1 DO *model_uv_arr[pol_i]+=*gal_model_uv[pol_i];*uv_mask_use
IF Min(Ptr_valid(gal_spectral_model_uv)) GT 0 THEN FOR pol_i=0,n_pol-1 DO FOR s_i=0,n_spectral-1 DO $
    *spectral_model_uv_arr[pol_i,s_i]+=*gal_spectral_model_uv[pol_i,s_i];*uv_mask_use
undefine_fhd,gal_model_uv,gal_spectral_model_uv

IF Keyword_Set(diffuse_filepath) THEN BEGIN
    IF file_test(diffuse_filepath) EQ 0 THEN diffuse_filepath=(file_search(diffuse_filepath+'*'))[0]
    print,"Reading diffuse model file: "+diffuse_filepath 
    diffuse_model_uv=fhd_diffuse_model(obs,jones,skymodel,spectral_model_arr=diffuse_spectral_model_uv,/uv_return,model_filepath=diffuse_filepath,_Extra=extra)
    IF Max(Ptr_valid(diffuse_model_uv)) EQ 0 THEN print,"Error reading or building diffuse model. Null pointer returned!"
ENDIF
IF Min(Ptr_valid(diffuse_model_uv)) GT 0 THEN FOR pol_i=0,n_pol-1 DO *model_uv_arr[pol_i]+=*diffuse_model_uv[pol_i];*uv_mask_use
IF Min(Ptr_valid(diffuse_spectral_model_uv)) GT 0 THEN FOR pol_i=0,n_pol-1 DO FOR s_i=0,n_spectral-1 DO $
    *spectral_model_uv_arr[pol_i,s_i]+=*diffuse_spectral_model_uv[pol_i,s_i];*uv_mask_use
undefine_fhd,diffuse_model_uv,diffuse_spectral_model_uv

vis_arr=Ptrarr(n_pol)

valid_test=fltarr(n_pol)
FOR pol_i=0,n_pol-1 DO valid_test[pol_i]=Total(Abs(*model_uv_arr[pol_i])) ; if the model only contains unpolarized sources but n_pol is set > 2, then this test will fail. Set n_pol=2.
IF min(valid_test) EQ 0 THEN BEGIN
    error=1
    print,"ERROR: Invalid model."
    timing=Systime(1)-t0
    RETURN,vis_arr
ENDIF

t_degrid=Fltarr(n_pol)
FOR pol_i=0,n_pol-1 DO BEGIN
    vis_arr[pol_i]=visibility_degrid(*model_uv_arr[pol_i],vis_weight_ptr[pol_i],obs,psf,params,silent=silent,$
        timing=t_degrid0,polarization=pol_i,fill_model_visibilities=fill_model_visibilities,$
        vis_input_ptr=vis_model_ptr[pol_i],spectral_model_uv_arr=spectral_model_uv_arr[pol_i,*], _Extra=extra)
    t_degrid[pol_i]=t_degrid0
ENDFOR
IF ~Keyword_Set(silent) THEN print,"Degridding timing: ",strn(t_degrid)

if keyword_set(model_delay_filter) then begin  
  ;params = replicate(params_original,1)
  vis_delay_filter, vis_arr, params, obs
  obs = replicate(obs_original,1)
endif
timing=Systime(1)-t0

RETURN,vis_arr
END
