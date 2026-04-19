module sfincs_openacc
   !
   use sfincs_data
   use sfincs_src_structures
   use sfincs_discharges,      only: qtsrc, nmindsrc
   use sfincs_rule_expression, only: rule_opcode, rule_atom, rule_cmp, rule_threshold, &
                                     rule_start, rule_length
   use sfincs_urban_drainage,  only: urban_drainage_zone_indices, urban_drainage_outfall_index, &
                                     urban_drainage_qmax, urban_drainage_backflow_coef, &
                                     urban_drainage_q_outfall, urban_drainage_cumulative_volume, &
                                     urb_zone_h_threshold, urb_zone_check_valve, &
                                     urb_zone_dh_design_min
   !
   implicit none
   !
contains
   !
   subroutine initialize_openacc()
   !
   !call acc_init( acc_device_nvidia )
   !
   ! Copy arrays to GPU memory
   ! 
   !$acc enter data, copyin( kcs, kfuv, kcuv, zs, zs0, zsderv, q, q0, uv, uv0, zb, zbuv, zbuvmx, zsmax, maxzsm, qmax, vmax, twet, zsm, z_volume, &
   !$acc               z_flags_iref, uv_flags_iref, uv_flags_type, uv_flags_dir, mask_adv, &
   !$acc               index_kcuv2, nmikcuv2, nmbkcuv2, ibkcuv2, zsb, zsb0, ibuvdir, uvmean, &
   !$acc               subgrid_uv_zmin, subgrid_uv_zmax, subgrid_uv_havg, subgrid_uv_nrep, subgrid_uv_pwet, &
   !$acc               subgrid_uv_havg_zmax, subgrid_uv_nrep_zmax, subgrid_uv_fnfit, subgrid_uv_navg_w, &
   !$acc               subgrid_z_zmin,  subgrid_z_zmax, subgrid_z_dep, subgrid_z_volmax, &
   !$acc               z_index_uv_md, z_index_uv_nd, z_index_uv_mu, z_index_uv_nu, &
   !$acc               uv_index_z_nm, uv_index_z_nmu, uv_index_u_nmd, uv_index_u_nmu, uv_index_u_ndm, uv_index_u_num, &
   !$acc               uv_index_v_ndm, uv_index_v_ndmu, uv_index_v_nm, uv_index_v_nmu, &
   !$acc               qsrc, qtsrc, q_src_struc, nmindsrc, src_struc_nm_in, src_struc_nm_out, src_struc_type, &
   !$acc               src_struc_direction, &
   !$acc               src_struc_nm_obs_1, src_struc_nm_obs_2, &
   !$acc               src_struc_q, src_struc_flow_coef, &
   !$acc               src_struc_width, src_struc_sill_elevation, src_struc_mannings_n, &
   !$acc               src_struc_opening_duration, src_struc_closing_duration, &
   !$acc               src_struc_height, &
   !$acc               src_struc_invert_1, src_struc_invert_2, &
   !$acc               src_struc_submergence_ratio, &
   !$acc               src_struc_distance, src_struc_status, src_struc_fraction_open, src_struc_t_state, &
   !$acc               src_struc_rule_open, src_struc_rule_close, &
   !$acc               rule_opcode, rule_atom, rule_cmp, rule_threshold, rule_start, rule_length, &
   !$acc               z_index_wavemaker, wavemaker_uvmean, wavemaker_nmd, wavemaker_nmu, wavemaker_ndm, wavemaker_num, &
   !$acc               structure_uv_index, structure_parameters, structure_type, structure_length, &
   !$acc               fwuv, &
   !$acc               tauwu, tauwv, tauwu0, tauwv0, tauwu1, tauwv1, &
   !$acc               windu, windv, windu0, windv0, windu1, windv1, windmax, &
   !$acc               patm, patm0, patm1, patmb, nmindbnd, &
   !$acc               prcp, prcp0, prcp1, cumprcp, netprcp, prcp, qext, &
   !$acc               dxminv, dxrinv, dyrinv, dxm2inv, dxr2inv, dyr2inv, dxrinvc, dyrinvc, dxm, dxrm, dyrm, cell_area_m2, cell_area, &
   !$acc               gn2uv, fcorio2d, storage_volume, nuvisc, &
   !$acc               cuv_index_uv, cuv_index_uv1, cuv_index_uv2, &
   !$acc               x73, &
   !$acc               gnapp2, &
   !$acc               timestep_analysis_required_timestep, timestep_analysis_average_required_timestep, timestep_analysis_times_wet, timestep_analysis_times_limiting, &
   !$acc               qinffield, qinfmap, cuminf, scs_rain, scs_Se, scs_P1, scs_F1, scs_S1, rain_T1, &
   !$acc               ksfield, GA_head, GA_sigma, GA_sigma_max, GA_F, GA_Lu, inf_kr, horton_kd, horton_fc, horton_f0, &
   !$acc               qdrain_rate, bucket_volume, bucket_capacity, bucket_k, bucket_drain_rate, bucket_loss, bucket_runoff, &
   !$acc               urban_drainage_zone_indices, urban_drainage_outfall_index, urban_drainage_qmax, urban_drainage_backflow_coef, &
   !$acc               urban_drainage_q_outfall, urban_drainage_cumulative_volume, &
   !$acc               urb_zone_h_threshold, urb_zone_check_valve, urb_zone_dh_design_min )
   !
   end subroutine
   !
   subroutine finalize_openacc()
   !
   !$acc exit data delete( kcs, kfuv, kcuv, zs, zs0, zsderv, q, q0, uv, uv0, zb, zbuv, zbuvmx, zsmax, maxzsm, qmax, vmax, twet, zsm, z_volume, &
   !$acc               z_flags_iref, uv_flags_iref, uv_flags_type, uv_flags_dir, mask_adv, &
   !$acc               index_kcuv2, nmikcuv2, nmbkcuv2, ibkcuv2, zsb, zsb0, ibuvdir, uvmean, &
   !$acc               subgrid_uv_zmin, subgrid_uv_zmax, subgrid_uv_havg, subgrid_uv_nrep, subgrid_uv_pwet, &
   !$acc               subgrid_uv_havg_zmax, subgrid_uv_nrep_zmax, subgrid_uv_fnfit, subgrid_uv_navg_w, &
   !$acc               subgrid_z_zmin,  subgrid_z_zmax, subgrid_z_dep, subgrid_z_volmax, &
   !$acc               z_index_uv_md, z_index_uv_nd, z_index_uv_mu, z_index_uv_nu, &
   !$acc               uv_index_z_nm, uv_index_z_nmu, uv_index_u_nmd, uv_index_u_nmu, uv_index_u_ndm, uv_index_u_num, &
   !$acc               uv_index_v_ndm, uv_index_v_ndmu, uv_index_v_nm, uv_index_v_nmu, &
   !$acc               qsrc, qtsrc, q_src_struc, nmindsrc, src_struc_nm_in, src_struc_nm_out, src_struc_type, &
   !$acc               src_struc_direction, &
   !$acc               src_struc_nm_obs_1, src_struc_nm_obs_2, &
   !$acc               src_struc_q, src_struc_flow_coef, &
   !$acc               src_struc_width, src_struc_sill_elevation, src_struc_mannings_n, &
   !$acc               src_struc_opening_duration, src_struc_closing_duration, &
   !$acc               src_struc_height, &
   !$acc               src_struc_invert_1, src_struc_invert_2, &
   !$acc               src_struc_submergence_ratio, &
   !$acc               src_struc_distance, src_struc_status, src_struc_fraction_open, src_struc_t_state, &
   !$acc               src_struc_rule_open, src_struc_rule_close, &
   !$acc               rule_opcode, rule_atom, rule_cmp, rule_threshold, rule_start, rule_length, &
   !$acc               z_index_wavemaker, wavemaker_uvmean, wavemaker_nmd, wavemaker_nmu, wavemaker_ndm, wavemaker_num, &
   !$acc               structure_uv_index, structure_parameters, structure_type, structure_length, &
   !$acc               fwuv, &
   !$acc               tauwu, tauwv, tauwu0, tauwv0, tauwu1, tauwv1, &
   !$acc               windu, windv, windu0, windv0, windu1, windv1, windmax, & 
   !$acc               patm, patm0, patm1, patmb, nmindbnd, &
   !$acc               prcp, prcp0, prcp1, cumprcp, qext, &
   !$acc               dxminv, dxrinv, dyrinv, dxm2inv, dxr2inv, dyr2inv, dxrinvc, dxm, dxrm, dyrm, cell_area_m2, cell_area, &
   !$acc               gn2uv, fcorio2d, storage_volume, nuvisc, &
   !$acc               cuv_index_uv, cuv_index_uv1, cuv_index_uv2, &
   !$acc               x73, &
   !$acc               gnapp2, &
   !$acc               timestep_analysis_required_timestep, timestep_analysis_average_required_timestep, timestep_analysis_times_wet, timestep_analysis_times_limiting, &   
   !$acc               qinffield, qinfmap, cuminf, scs_rain, scs_Se, scs_P1, scs_F1, scs_S1, rain_T1, &
   !$acc               ksfield, GA_head, GA_sigma, GA_sigma_max, GA_F, GA_Lu, inf_kr, horton_kd, horton_fc, horton_f0, &
   !$acc               qdrain_rate, bucket_volume, bucket_capacity, bucket_k, bucket_drain_rate, bucket_loss, bucket_runoff, &
   !$acc               urban_drainage_zone_indices, urban_drainage_outfall_index, urban_drainage_qmax, urban_drainage_backflow_coef, &
   !$acc               urban_drainage_q_outfall, urban_drainage_cumulative_volume, &
   !$acc               urb_zone_h_threshold, urb_zone_check_valve, urb_zone_dh_design_min )
   !
   end
   !
end module
