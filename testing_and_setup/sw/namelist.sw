&sw_model
    config_test_case = 5
    config_time_integration = 'LTS3'
    config_dt = 200
    config_calendar_type = 'gregorian_noleap'
    config_start_time = '0001-01-01_00:00:00'
    config_stop_time = 'none'
    config_run_duration = '00:05:00'
    config_stats_interval = 100
    config_h_ScaleWithMesh = false
    config_h_mom_eddy_visc2 = 0.0
    config_h_mom_eddy_visc4 = 0.0
    config_h_tracer_eddy_diff2 = 0.0
    config_h_tracer_eddy_diff4 = 0.0
    config_thickness_adv_order = 2
    config_tracer_adv_order = 2
    config_positive_definite = false
    config_monotonic = false
    config_wind_stress = false
    config_bottom_drag = false
    config_apvm_upwinding = 0.0
    config_num_halos = 2
/
&io
    config_pio_num_iotasks = 0
    config_pio_stride = 1
/
&decomposition
    config_block_decomp_file_prefix = 'graph.info.part.'
    config_number_of_blocks = 12
    config_explicit_proc_decomp = false
    config_proc_decomp_file_prefix = 'block_proc_list.part.'
/
&restart
    config_do_restart = false
/
&local_time_stepping
   config_use_local_time_stepping = true
   config_dt_scaling_LTS = 25
/
