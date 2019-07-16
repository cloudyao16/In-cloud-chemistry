!4 Copyright (C) 2009-2013 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The process program.

!> Read NetCDF output files and process them.
program process

  use pmc_output
  use pmc_stats

  character(len=PMC_MAX_FILENAME_LEN), parameter :: prefix &
       = "scenarios/mono2/out/urban_plume_aq_chem"

  character(len=PMC_MAX_FILENAME_LEN) :: in_filename, out_filename
  type(bin_grid_t)   :: diam_grid, bc_grid, sc_grid, su_grid, ph_grid,entropy_grid, avg_bin_grid
  type(aero_data_t)  :: aero_data
  type(aero_state_t) :: aero_state, aero_state_averaged
  type(env_state_t)  :: env_state
  !type(gas_state_t) :: gas_state
  integer :: ncid, index, repeat, i_index, i_repeat, n_index, n_repeat
  logical, allocatable :: sr_neut(:), sr_acid(:)
  real(kind=dp) :: time, del_t, tot_num_conc, tot_mass_conc, totdry_mass_conc,&
                   tot_entropy,tot_bc_conc,tot_so4_conc,&
                   tot_oc_conc,tot_nh4_conc,tot_h2o_conc,rh,temp,spres,pre,density,&
                   tot_aO3_conc, tot_aSO2_conc, tot_aHCHO_conc, & ! sulfate oxidation
                   tot_HSO4m_conc, tot_HSO3m_conc, tot_SO3m_conc, tot_SO4m_conc, &
                   tot_SO3mm_conc, tot_SO5m_conc, tot_HSO5m_conc, tot_SO5O2Hm_conc,&
                   tot_SO5O2mm_conc, tot_Hp_conc, tot_FEpp_conc,& ! Sulfate all added 
                   tot_BRm_conc, tot_aHNO4_conc, & 
                   tot_no3_conc, tot_aNO3_conc, tot_aN2O5_conc, tot_aNO_conc,&
                   tot_NO2p_conc, tot_NO4m_conc, tot_NO2m_conc, tot_OHm_conc, tot_aH2O2_conc, &
                   tot_aHO2_conc, tot_aNO2_conc, tot_aNH3_conc, tot_aOH_conc,&
                   tot_aH2C2O4_conc, tot_C2O4m_conc, tot_HC2O4m_conc, &
                   tot_FEppp_conc, tot_FEOHpp_conc, tot_FEC2O4p_conc, tot_aCO2_conc, &
                   tot_all_conc, tot_aO2_conc
  real(kind=dp) :: tot_entropy_averaged
  character(len=PMC_UUID_LEN) :: uuid
  real(kind=dp), allocatable :: times(:), dry_diameters(:), wet_diameters(:),num_concs(:), &
       dry_masses(:), masses(:), bc_masses(:), bc_fracs(:),oc_masses(:), &
       so4_masses(:),su_fracs(:),no3_masses(:),nh4_masses(:),h2o_masses(:),&
       num_concs_averaged(:), dry_masses_averaged(:), masses_averaged(:), &
       entropies(:), entropies_averaged(:), crit_rhs(:), scs(:),&
       num_dist(:),bc_dist_dry(:), su_dist_dry_acid(:), su_dist_dry_neut(:),&
       oc_dist_dry(:),no3_dist_dry(:),&
       nh4_dist_dry_acid(:), nh4_dist_dry_neut(:), tot_dist_dry(:),&
       bc_dist_wet(:),su_dist_wet(:),oc_dist_wet(:),no3_dist_wet(:),nh4_dist_wet(:), tot_dist_wet(:),&
       diam_bc_dist(:,:), diam_sc_dist(:,:), diam_su_dist(:,:), diam_su_dist_wet(:,:),diam_ph_dist(:,:), &
       entropy_dist(:),wet_num_dist(:), diam_ph_dist_acid(:,:), diam_ph_dist_neut(:,:), &
       aO3_masses(:), aSO2_masses(:), aHCHO_masses(:), HSO4m_masses(:), HSO3m_masses(:), &
       SO3m_masses(:), SO4m_masses(:), SO3mm_masses(:), SO5m_masses(:), HSO5m_masses(:), &
       SO5O2Hm_masses(:), SO5O2mm_masses(:), Hp_masses(:), FEpp_masses(:),&
       BRm_masses(:), aHNO4_masses(:), aH2O2_masses(:),&
       aNO3_masses(:), aN2O5_masses(:), aNO_masses(:), NO2p_masses(:), NO4m_masses(:),&
       NO2m_masses(:), aHO2_masses(:), aNO2_masses(:), aNH3_masses(:), aOH_masses(:),&
       Hp_dist_dry(:),Hp_dist_wet(:),h2o_dist_dry(:),h2o_dist_wet(:),ph(:), ph_acid(:), ph_neut(:), &
       su_frac_acid(:), su_frac_neut(:), nh4_frac_acid(:), nh4_frac_neut(:), &
       aH2C2O4_masses(:), C2O4m_masses(:), HC2O4m_masses(:), FEppp_masses(:),&
       FEOHpp_masses(:), FEC2O4p_masses(:), lwc(:), neut(:), acid(:), aCO2_masses(:), all_masses(:), &
       aO2_masses(:), Hp_acid(:), Hp_neut(:), h2o_acid(:), h2o_neut(:), Na_masses(:) 
       !su_num_dist(:), bc_num_dist(:), oc_num_dist(:), no3_num_dist(:),&
       !su_num_concs(:),bc_num_concs(:), oc_num_concs(:), no3_num_concs(:)
  type(stats_1d_t) :: stats_num_dist, stats_wet_num_dist, stats_entropy_dist, stats_tot_num_conc, &
       stats_tot_mass_conc, stats_totdry_mass_conc,stats_tot_entropy, stats_tot_entropy_averaged, &
       stats_tot_bc_conc,stats_tot_so4_conc,stats_tot_no3_conc,stats_tot_nh4_conc,&
       stats_tot_oc_conc,stats_tot_h2o_conc,stats_rh,stats_temp,stats_spres,stats_pre,stats_den,&
       stats_bc_dist_dry,stats_su_dist_dry,stats_oc_dist_dry,stats_no3_dist_dry,stats_nh4_dist_dry_acid,&
       stats_nh4_dist_dry_neut, stats_su_dist_dry_acid, stats_su_dist_dry_neut, &
       stats_bc_dist_wet,stats_su_dist_wet,stats_oc_dist_wet,stats_no3_dist_wet,stats_nh4_dist_wet, &
       stats_tot_aO3_conc, stats_tot_aSO2_conc, stats_tot_aHCHO_conc, stats_tot_HSO4m_conc, & 
       stats_tot_HSO3m_conc, stats_tot_SO3m_conc, stats_tot_SO4m_conc, &
       stats_tot_SO3mm_conc, stats_tot_SO5m_conc, stats_tot_HSO5m_conc, &
       stats_tot_SO5O2Hm_conc, stats_tot_SO5O2mm_conc, stats_tot_Hp_conc,&! sulfate all added
       stats_tot_FEpp_conc, stats_tot_FEppp_conc, stats_tot_FEOHpp_conc, stats_tot_FEC2O4p_conc,&
       stats_tot_BRm_conc, stats_tot_aHNO4_conc, stats_tot_aNO3_conc,&
       stats_tot_aN2O5_conc, stats_tot_aNO_conc,& 
       stats_tot_NO2p_conc, stats_tot_NO4m_conc, stats_tot_NO2m_conc, &
       stats_tot_aH2O2_conc, stats_tot_aHO2_conc, stats_tot_aNO2_conc, &
       stats_tot_aNH3_conc, stats_tot_aOH_conc,&
       stats_Hp_dist_dry, stats_Hp_dist_wet, stats_h2o_dist_dry, stats_h2o_dist_wet, &
       stats_tot_aH2C2O4_conc, stats_tot_C2O4m_conc, stats_tot_HC2O4m_conc, &
       stats_tot_aCO2_conc, stats_tot_aO2_conc
type(stats_2d_t) :: stats_diam_bc_dist, stats_diam_sc_dist,stats_diam_su_dist,& 
stats_diam_ph_dist, stats_diam_su_dist_wet, stats_diam_ph_dist_acid,&
stats_diam_ph_dist_neut

  call pmc_mpi_init()

  call bin_grid_allocate(diam_grid)
  call bin_grid_allocate(bc_grid)
  call bin_grid_allocate(su_grid)
  call bin_grid_allocate(sc_grid)
  call bin_grid_allocate(ph_grid)
  call bin_grid_allocate(entropy_grid)
  call aero_data_allocate(aero_data)
  call aero_state_allocate(aero_state)
  call aero_state_allocate(aero_state_averaged)
  call bin_grid_allocate(avg_bin_grid)
  call env_state_allocate(env_state)
  !call gas_state_allocate(gas_state)
  call input_n_files(prefix, n_repeat, n_index)

  call bin_grid_make(diam_grid, BIN_GRID_TYPE_LOG, 180, 1d-9, 1d-3)
  call bin_grid_make(bc_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(su_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(ph_grid, BIN_GRID_TYPE_LINEAR, 20, 1d0, 1d1)
  call bin_grid_make(sc_grid, BIN_GRID_TYPE_LOG, 50, 1d-4, 1d0)
  call bin_grid_make(entropy_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(avg_bin_grid, BIN_GRID_TYPE_LOG, 1, 1d-30, 1d10)

  allocate(times(n_index))

  do i_index = 1,n_index
     do i_repeat = 1,n_repeat
        call make_filename(in_filename, prefix, ".nc", i_index, i_repeat)
        write(*,*) "Processing " // trim(in_filename)
        ! FIXME: add UUID check into input_state(), keyed off of index or
        ! time or something? Or init to "" and check if not this.
        call input_state(in_filename, index, time, del_t, repeat, &
             uuid, aero_data=aero_data, aero_state=aero_state, &
             env_state=env_state)!, gas_state=gas_state)

        times(i_index) = time

        dry_diameters = aero_state_dry_diameters(aero_state, aero_data)
        num_concs = aero_state_num_concs(aero_state)
        num_dist = bin_grid_histogram_1d(diam_grid, dry_diameters, num_concs)
        call stats_1d_add(stats_num_dist, num_dist)
        !write(*,*) size(num_dist) 
   
        ! Add the wet diameters
        wet_diameters = aero_state_diameters(aero_state)
        wet_num_dist  = bin_grid_histogram_1d(diam_grid, wet_diameters,num_concs)
        call stats_1d_add(stats_wet_num_dist, wet_num_dist)
        
       ! Add diagnose for wet diam
        max_dp = 0.0
        do i=1,size(wet_diameters)
            if (wet_diameters(i).gt.max_dp) then
                max_dp = wet_diameters(i)
            endif
        enddo

        ! add the source 
        neut  = aero_state_src(aero_state, aero_data, include=(/"init_neut"/)) 
        acid  = aero_state_src(aero_state, aero_data, include=(/"init_acid"/)) 
        sr_neut =  neut .gt. 0
        sr_acid =  acid .gt. 0 
        
        !do i = 1, size(neut)
        !   write(*, 1000) neut(i)
        !enddo
   !1000 format (F3.1)
        !write(*,*) "time ",times(i_index),"Max diameter ",max_dp
        tot_num_conc = sum(num_concs)
        call stats_1d_add_entry(stats_tot_num_conc, tot_num_conc, i_index)

        masses = aero_state_masses(aero_state, aero_data)
        tot_mass_conc = sum(masses * num_concs)
        call stats_1d_add_entry(stats_tot_mass_conc, tot_mass_conc, i_index)

        dry_masses = aero_state_masses(aero_state, aero_data, &
             exclude=(/"H2O"/))
        totdry_mass_conc = sum(dry_masses * num_concs)
        call stats_1d_add_entry(stats_totdry_mass_conc, totdry_mass_conc, i_index)
 
        h2o_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"H2O"/))
        tot_h2o_conc = sum(h2o_masses*num_concs)
        call stats_1d_add_entry(stats_tot_h2o_conc, tot_h2o_conc, i_index)
        h2o_dist_dry = bin_grid_histogram_1d(diam_grid, dry_diameters,h2o_masses*num_concs)
        h2o_dist_wet = bin_grid_histogram_1d(diam_grid, wet_diameters,h2o_masses*num_concs)
        call stats_1d_add(stats_h2o_dist_dry, h2o_dist_dry)
        call stats_1d_add(stats_h2o_dist_wet, h2o_dist_wet)
        
        bc_masses   = aero_state_masses(aero_state, aero_data, &
             include=(/"BC"/))
        tot_bc_conc = sum(bc_masses*num_concs)
        call stats_1d_add_entry(stats_tot_bc_conc, tot_bc_conc, i_index)
        !write(*,*) size(bc_masses) 
        bc_dist_dry = bin_grid_histogram_1d(diam_grid, dry_diameters,bc_masses*num_concs)
        bc_dist_wet = bin_grid_histogram_1d(diam_grid, wet_diameters,bc_masses*num_concs)
        call stats_1d_add(stats_bc_dist_dry, bc_dist_dry)
        call stats_1d_add(stats_bc_dist_wet, bc_dist_wet)

        oc_masses = aero_state_masses(aero_state, aero_data, &
             include = (/"OC"/))
        tot_oc_conc = sum(oc_masses*num_concs)
        call stats_1d_add_entry(stats_tot_oc_conc, tot_oc_conc, i_index)
        oc_dist_dry = bin_grid_histogram_1d(diam_grid, dry_diameters, oc_masses*num_concs)
        oc_dist_wet = bin_grid_histogram_1d(diam_grid, wet_diameters, oc_masses*num_concs)
        call stats_1d_add(stats_oc_dist_dry, oc_dist_dry)
        call stats_1d_add(stats_oc_dist_wet, oc_dist_wet)        

        so4_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"SO4  ", "HSO4m"/))
        tot_so4_conc     = sum(so4_masses*num_concs)
        
        call stats_1d_add_entry(stats_tot_so4_conc, tot_so4_conc, i_index)
        su_dist_dry_acid = bin_grid_histogram_1d(diam_grid, pack(dry_diameters,sr_acid), &
        pack(so4_masses*num_concs, sr_acid))
        call stats_1d_add(stats_su_dist_dry_acid, su_dist_dry_acid)
   
        su_dist_dry_neut = bin_grid_histogram_1d(diam_grid, pack(dry_diameters,sr_neut), &
        pack(so4_masses*num_concs, sr_neut))
        call stats_1d_add(stats_su_dist_dry_neut, su_dist_dry_neut)
        
        su_dist_wet = bin_grid_histogram_1d(diam_grid, wet_diameters, so4_masses*num_concs)
        call stats_1d_add(stats_su_dist_wet, su_dist_wet)

        nh4_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"NH4"/))
        tot_nh4_conc = sum(nh4_masses*num_concs)
        call stats_1d_add_entry(stats_tot_nh4_conc, tot_nh4_conc, i_index)
        nh4_dist_dry_acid = bin_grid_histogram_1d(diam_grid, &
        pack(dry_diameters,sr_acid), pack(nh4_masses*num_concs,sr_acid))
        call stats_1d_add(stats_nh4_dist_dry_acid, nh4_dist_dry_acid)

        nh4_dist_dry_neut = bin_grid_histogram_1d(diam_grid, &
        pack(dry_diameters,sr_neut), pack(nh4_masses*num_concs,sr_neut))
        call stats_1d_add(stats_nh4_dist_dry_neut, nh4_dist_dry_neut)

        nh4_dist_wet = bin_grid_histogram_1d(diam_grid, wet_diameters, nh4_masses*num_concs )
        call stats_1d_add(stats_nh4_dist_wet, nh4_dist_wet)

        Hp_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"Hp"/))
        tot_Hp_conc  = sum(Hp_masses*num_concs)
        call stats_1d_add_entry(stats_tot_Hp_conc, tot_Hp_conc, i_index)
        Hp_dist_dry = bin_grid_histogram_1d(diam_grid, dry_diameters,Hp_masses*num_concs)
        Hp_dist_wet = bin_grid_histogram_1d(diam_grid, wet_diameters,Hp_masses*num_concs)
        call stats_1d_add(stats_Hp_dist_dry, Hp_dist_dry)
        call stats_1d_add(stats_Hp_dist_wet, Hp_dist_wet)
        
        ph  = -log10(Hp_masses*1000/h2o_masses)
        lwc =  h2o_masses*num_concs * 1000 ! kg/m3 to g/m3
        diam_ph_dist_acid = bin_grid_histogram_2d(diam_grid, pack(dry_diameters,sr_acid), & 
                            ph_grid, pack(ph, sr_acid), pack(num_concs,sr_acid))
        diam_ph_dist_neut = bin_grid_histogram_2d(diam_grid, pack(dry_diameters,sr_neut), & 
                            ph_grid, pack(ph, sr_neut), pack(num_concs,sr_neut))
        call stats_2d_add(stats_diam_ph_dist_acid, diam_ph_dist_acid)
        call stats_2d_add(stats_diam_ph_dist_neut, diam_ph_dist_neut)
   
        Na_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"Na"/))
 
              ph_acid = pack(ph, sr_acid) 
        nh4_frac_acid = pack(nh4_masses/dry_masses*100, sr_acid) 
        su_frac_acid  = pack(so4_masses/dry_masses*100, sr_acid) 

          ph_neut     = pack(ph, sr_neut)                       
        nh4_frac_neut = pack(nh4_masses/dry_masses*100, sr_neut)
        su_frac_neut  = pack(so4_masses/dry_masses*100, sr_neut)

           Hp_neut    = pack(Hp_masses, sr_neut)
          h2o_neut    = pack(h2o_masses, sr_neut)
           Hp_acid    = pack(Hp_masses, sr_acid)
          h2o_acid    = pack(h2o_masses, sr_acid)
 
        
        do i=1,size(ph_acid)   
           write(*,*) ph_acid(i), so4_masses(i)
        enddo

         !do i=1,size(ph_neut)
         !   write(*,*) ph_neut(i), so4_masses(i) 
         !enddo

        !do i=1,size(wet_diameters)
        !    write(*,*) ph(i), so4_masses(i), so4_masses(i)/dry_masses(i)*100
        !enddo

        !all_masses    = aero_state_masses(aero_state, aero_data, &
        !  include=(/"SO4    ","NO3    ","Cl     ","NH4    ","MSA    ","ARO1   ","ARO2   ","ALK1   ",&
        !            "OLE1   ","API1   ","API2   ","LIM1   ","LIM2   ","CO3    ","Na     ","Ca     ",&
        !            "OIN    ","OC     ","BC     ","aCO2   ","aO3    ","aHO2   ","aOH    ", &
                    !"HO2m   ","aH2O2  ","aqNO2  ","aHONO  ","aNO3   ","aN2O5  ","aNH3   ","aHCHO  ",&
                    !"aORA1  ","aSO2   ","aOP1   ","aORA2  ","aMO2   ","aETHPX ","aETOH  ","aCH3OH ",&
        !            "aALD   ","aBR2   ","aCL2   ","aHNO4  ","aACO3  ","aGLY   ","aO2    ","aCLNO2 "/))
                    !"aBRNO2 ","aBRCL  ","aNO    ","FEOHpp ","FEpp   ","OHm    ","C2O4mm ","CO2m   ",&
                    !"FEppp  ","CUp    ","CUpp   ","O2m    ","Hp     ","O3m    ","aHO3   ","HSO3m  ",&
                    !"SO3m   ","FEOpp  ","CLOHm  ","NO2p   ","SO4m   ","NO4m   ","NO2m   ","HSO4m  "/))
                    !"BRm    ","HMSm   ","CHOSO3m","aO2CHO ","SO3mm  ","SO5m   ","HSO5m  ","SO5O2Hm",&
                    !"SO5O2mm","aCH2OH ","aCH2OH2","aCHOH2 ","aCH3CO ","aCO2H  ","HCOOm  ","MCOOm  ",&
                    !"CH2COOm","aCH3O  ","HC2O4m ","C2O4m  ","aH2C2O4","CL2m   ","aHOCL  ","aBR    "/))
                    !"BR2m   ","aHOBR  ","BROHm  ","aCL    ","FEC2O4p","NAp_C  "/))
        all_masses    = aero_state_masses(aero_state, aero_data, &
          include=(/"SO4  ", "HSO3m", "aCO2 ", "aO2  ", "HSO4m", "NH4  "/))

        tot_all_conc     = sum(all_masses*num_concs)
        !write(*,*) tot_all_conc/totdry_mass_conc
        !%%%%%%%%%%%%SO4 RELATED SPECIES%%%%%%%%
        aO3_masses    = aero_state_masses(aero_state, aero_data, &
               include=(/"aO3"/))
        tot_aO3_conc     = sum(aO3_masses*num_concs)
        call stats_1d_add_entry(stats_tot_aO3_conc, tot_aO3_conc, i_index)

        aH2O2_masses    = aero_state_masses(aero_state, aero_data, &
               include=(/"aH2O2"/))
        tot_aH2O2_conc     = sum(aH2O2_masses*num_concs)
        call stats_1d_add_entry(stats_tot_aH2O2_conc, tot_aH2O2_conc, i_index)
        
        aNO2_masses    = aero_state_masses(aero_state, aero_data, &
               include=(/"aqNO2"/))
        tot_aNO2_conc     = sum(aNO2_masses*num_concs)
        call stats_1d_add_entry(stats_tot_aNO2_conc, tot_aNO2_conc, i_index)
        
        aNH3_masses    = aero_state_masses(aero_state, aero_data, &
               include=(/"aNH3"/))
        tot_aNH3_conc     = sum(aNH3_masses*num_concs)
        call stats_1d_add_entry(stats_tot_aNH3_conc, tot_aNH3_conc, i_index)

        aOH_masses    = aero_state_masses(aero_state, aero_data, &
               include=(/"aOH"/))
        tot_aOH_conc     = sum(aOH_masses*num_concs)
        call stats_1d_add_entry(stats_tot_aOH_conc, tot_aOH_conc, i_index)

        aHO2_masses    = aero_state_masses(aero_state, aero_data, &
               include=(/"aHO2"/))
        tot_aHO2_conc  = sum(aHO2_masses*num_concs)
        call stats_1d_add_entry(stats_tot_aHO2_conc, tot_aHO2_conc, i_index)

        FEpp_masses    = aero_state_masses(aero_state, aero_data, &
               include=(/"FEpp"/))
        tot_FEpp_conc  = sum(FEpp_masses*num_concs)
        call stats_1d_add_entry(stats_tot_FEpp_conc, tot_FEpp_conc, i_index)

        FEppp_masses    = aero_state_masses(aero_state, aero_data, &
               include=(/"FEppp"/))
        tot_FEppp_conc  = sum(FEppp_masses*num_concs)
        call stats_1d_add_entry(stats_tot_FEppp_conc, tot_FEppp_conc, i_index)

        FEOHpp_masses    = aero_state_masses(aero_state, aero_data, &
               include=(/"FEOHpp"/))
        tot_FEOHpp_conc  = sum(FEOHpp_masses*num_concs)
        call stats_1d_add_entry(stats_tot_FEOHpp_conc, tot_FEOHpp_conc, i_index)

        FEC2O4p_masses    = aero_state_masses(aero_state, aero_data, &
               include=(/"FEC2O4p"/))
        tot_FEC2O4p_conc  = sum(FEC2O4p_masses*num_concs)
        call stats_1d_add_entry(stats_tot_FEC2O4p_conc, tot_FEC2O4p_conc, i_index)
     
        BRm_masses    = aero_state_masses(aero_state, aero_data, &
               include=(/"BRm"/))
        tot_BRm_conc  = sum(BRm_masses*num_concs)
        call stats_1d_add_entry(stats_tot_BRm_conc, tot_BRm_conc, i_index)

        aHNO4_masses    = aero_state_masses(aero_state, aero_data, &
               include=(/"aHNO4"/))
        tot_aHNO4_conc  = sum(aHNO4_masses*num_concs) 
        call stats_1d_add_entry(stats_tot_aHNO4_conc, tot_aHNO4_conc, i_index)

        aSO2_masses   = aero_state_masses(aero_state, aero_data, &
               include=(/"aSO2"/))
        tot_aSO2_conc  = sum(aSO2_masses*num_concs)
        call stats_1d_add_entry(stats_tot_aSO2_conc, tot_aSO2_conc, i_index)

        aHCHO_masses  = aero_state_masses(aero_state, aero_data, &
               include=(/"aHCHO"/))
        tot_aHCHO_conc= sum(aHCHO_masses*num_concs)
        call stats_1d_add_entry(stats_tot_aHCHO_conc, tot_aHCHO_conc, i_index)

        HSO4m_masses  = aero_state_masses(aero_state, aero_data, &
               include=(/"HSO4m"/))
        tot_HSO4m_conc= sum(HSO4m_masses*num_concs)
        call stats_1d_add_entry(stats_tot_HSO4m_conc, tot_HSO4m_conc, i_index)
       
        HSO3m_masses  = aero_state_masses(aero_state, aero_data, &
               include=(/"HSO3m"/))
        tot_HSO3m_conc    = sum(HSO3m_masses*num_concs)
        call stats_1d_add_entry(stats_tot_HSO3m_conc, tot_HSO3m_conc, i_index)

        SO3m_masses   = aero_state_masses(aero_state, aero_data, &
               include=(/"SO3m"/))
        tot_SO3m_conc    = sum(SO3m_masses*num_concs)
        call stats_1d_add_entry(stats_tot_SO3m_conc, tot_SO3m_conc, i_index)

        SO4m_masses   = aero_state_masses(aero_state, aero_data, &
               include=(/"SO4m"/))
        tot_SO4m_conc= sum(SO4m_masses*num_concs)
        call stats_1d_add_entry(stats_tot_SO4m_conc, tot_SO4m_conc, i_index)

        SO3mm_masses  = aero_state_masses(aero_state, aero_data, &
               include=(/"SO3mm"/))
        tot_SO3mm_conc= sum(SO3mm_masses*num_concs)
        call stats_1d_add_entry(stats_tot_SO3mm_conc, tot_SO3mm_conc, i_index)

        SO5m_masses  = aero_state_masses(aero_state, aero_data, &
               include=(/"SO5m"/))
        tot_SO5m_conc= sum(SO5m_masses*num_concs)
        call stats_1d_add_entry(stats_tot_SO5m_conc, tot_SO5m_conc, i_index)

        HSO5m_masses  = aero_state_masses(aero_state, aero_data, &
               include=(/"HSO5m"/))
        tot_HSO5m_conc= sum(HSO5m_masses*num_concs)
        call stats_1d_add_entry(stats_tot_HSO5m_conc, tot_HSO5m_conc, i_index)

        SO5O2Hm_masses= aero_state_masses(aero_state, aero_data, &
               include=(/"SO5O2Hm"/))
        tot_SO5O2Hm_conc= sum(SO5O2Hm_masses*num_concs)
        call stats_1d_add_entry(stats_tot_SO5O2Hm_conc, tot_SO5O2Hm_conc, i_index)

        SO5O2mm_masses= aero_state_masses(aero_state, aero_data, &
               include=(/"SO5O2mm"/))
        tot_SO5O2mm_conc = sum(SO5O2mm_masses*num_concs)
        call stats_1d_add_entry(stats_tot_SO5O2mm_conc, tot_SO5O2mm_conc, i_index)
        !%%%%%%%%%%%%SO4 RELATED SPECIES%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%Nitric acid RELATED SPECIES%%%%%%%%
        aNO3_masses= aero_state_masses(aero_state, aero_data, &
               include=(/"aNO3"/))
        tot_aNO3_conc = sum(aNO3_masses*num_concs)
        call stats_1d_add_entry(stats_tot_aNO3_conc, tot_aNO3_conc, i_index) 

        aN2O5_masses= aero_state_masses(aero_state, aero_data, &
               include=(/"aN2O5"/))
        tot_aN2O5_conc = sum(aN2O5_masses*num_concs)
        call stats_1d_add_entry(stats_tot_aN2O5_conc, tot_aN2O5_conc, i_index)

        aNO_masses= aero_state_masses(aero_state, aero_data, &
               include=(/"aNO"/))
        tot_aNO_conc = sum(aNO_masses*num_concs)
        call stats_1d_add_entry(stats_tot_aNO_conc, tot_aNO_conc, i_index)
        
        NO2p_masses= aero_state_masses(aero_state, aero_data, &
               include=(/"NO2p"/))
        tot_NO2p_conc = sum(NO2p_masses*num_concs)
        call stats_1d_add_entry(stats_tot_NO2p_conc, tot_NO2p_conc, i_index)
        
        NO4m_masses= aero_state_masses(aero_state, aero_data, &
               include=(/"NO4m"/))
        tot_NO4m_conc = sum(NO4m_masses*num_concs)
        call stats_1d_add_entry(stats_tot_NO4m_conc, tot_NO4m_conc, i_index)

        NO2m_masses= aero_state_masses(aero_state, aero_data, &
               include=(/"NO2m"/))
        tot_NO2m_conc = sum(NO2m_masses*num_concs)
        call stats_1d_add_entry(stats_tot_NO2m_conc, tot_NO2m_conc, i_index)
        !%%%%%%%%%%%%Nitric RELATED SPECIES%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         
        !%%%%%%%%For test the negative values%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        aH2C2O4_masses= aero_state_masses(aero_state, aero_data, &
               include=(/"aH2C2O4"/))
        tot_aH2C2O4_conc = sum(aH2C2O4_masses*num_concs)
        call stats_1d_add_entry(stats_tot_aH2C2O4_conc, tot_aH2C2O4_conc, i_index)

        C2O4m_masses= aero_state_masses(aero_state, aero_data, &
               include=(/"C2O4m"/))
        tot_C2O4m_conc = sum(C2O4m_masses*num_concs)
        call stats_1d_add_entry(stats_tot_C2O4m_conc, tot_C2O4m_conc, i_index)

        HC2O4m_masses= aero_state_masses(aero_state, aero_data, &
               include=(/"HC2O4m"/))
        tot_HC2O4m_conc = sum(HC2O4m_masses*num_concs)
        call stats_1d_add_entry(stats_tot_HC2O4m_conc, tot_HC2O4m_conc, i_index)
       !%%%%%%%%%%%%%%%NEGATIVE VALUES TEST%%%%%%%%%%%%%%
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        !%%%%%%%%aCO2 and aO2 species%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        aCO2_masses= aero_state_masses(aero_state, aero_data, &
               include=(/"aCO2"/))
        tot_aCO2_conc = sum(aCO2_masses*num_concs)
        call stats_1d_add_entry(stats_tot_aCO2_conc, tot_aCO2_conc, i_index)

        aO2_masses= aero_state_masses(aero_state, aero_data, &
               include=(/"aO2"/))
        tot_aO2_conc = sum(aO2_masses*num_concs)
        call stats_1d_add_entry(stats_tot_aO2_conc, tot_aO2_conc, i_index)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        no3_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"NO3"/))
        tot_no3_conc = sum(no3_masses*num_concs)
        call stats_1d_add_entry(stats_tot_no3_conc, tot_no3_conc, i_index)
        no3_dist_dry = bin_grid_histogram_1d(diam_grid, dry_diameters, no3_masses*num_concs)
        no3_dist_wet = bin_grid_histogram_1d(diam_grid, wet_diameters, no3_masses*num_concs)
        call stats_1d_add(stats_no3_dist_dry, no3_dist_dry)
        call stats_1d_add(stats_no3_dist_wet, no3_dist_wet)

        bc_fracs = bc_masses / dry_masses
        diam_bc_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             bc_grid, bc_fracs, num_concs)
        call stats_2d_add(stats_diam_bc_dist, diam_bc_dist)

        su_fracs = so4_masses / dry_masses
        diam_su_dist     = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             su_grid, su_fracs, num_concs)
        diam_su_dist_wet = bin_grid_histogram_2d(diam_grid, wet_diameters, &
             su_grid, su_fracs, num_concs)
        call stats_2d_add(stats_diam_su_dist_wet, diam_su_dist_wet)
        call stats_2d_add(stats_diam_su_dist, diam_su_dist)

        crit_rhs = aero_state_crit_rel_humids(aero_state, aero_data, &
             env_state)
        scs = crit_rhs - 1d0
        diam_sc_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             sc_grid, scs, num_concs)
        call stats_2d_add(stats_diam_sc_dist, diam_sc_dist)

        entropies = aero_state_mass_entropies(aero_state, aero_data) !, &
             !exclude=["H2O"]) !, group=["BC"])
        entropy_dist = bin_grid_histogram_1d(entropy_grid, entropies, &
             num_concs)
        call stats_1d_add(stats_entropy_dist, entropy_dist)

        tot_entropy = sum(entropies * masses * num_concs) &
             / sum(masses * num_concs)
        call stats_1d_add_entry(stats_tot_entropy, tot_entropy, i_index)
      
        ! Add RH, temperature,density & pressure
        rh   = env_state%rel_humid
        call stats_1d_add_entry(stats_rh,rh,i_index)

        temp = env_state%temp
        call stats_1d_add_entry(stats_temp,temp,i_index)

        spres= env_state_sat_vapor_pressure(env_state)
        call stats_1d_add_entry(stats_spres,spres,i_index)
       
        pre  = env_state%pressure
        call stats_1d_add_entry(stats_pre,pre,i_index)

        density = env_state_air_den(env_state)
        call stats_1d_add_entry(stats_den, density, i_index)
         
        ! Add the gas variables
        
        call aero_state_copy(aero_state, aero_state_averaged)
        call aero_state_bin_average_comp(aero_state_averaged, avg_bin_grid, &
             aero_data)
        num_concs_averaged = aero_state_num_concs(aero_state_averaged)
        masses_averaged = aero_state_masses(aero_state_averaged, aero_data)
        dry_masses_averaged = aero_state_masses(aero_state_averaged, &
             aero_data, exclude=(/"H2O"/))
        entropies_averaged = aero_state_mass_entropies(aero_state_averaged, &
             aero_data) !, exclude=["H2O"]) !, group=["BC"])
        tot_entropy_averaged &
             = sum(entropies_averaged * masses_averaged &
             * num_concs_averaged) &
             / sum(masses_averaged * num_concs_averaged)
        call stats_1d_add_entry(stats_tot_entropy_averaged, &
             tot_entropy_averaged, i_index)
        
         !do i=1,size(dry_diameters)
         !      write(*,*)  i_index, dry_diameters(i) , wet_diameters(i), crit_rhs(i), ph(i),& 
         !                  lwc(i), so4_masses(i), dry_masses(i), masses(i)
         !enddo
       
     end do

     call make_filename(out_filename, prefix, "_process.nc", index)
     call pmc_nc_open_write(out_filename, ncid)
     call pmc_nc_write_info(ncid, uuid, "8_urban_plume_aq_chem process")
     call bin_grid_output_netcdf(diam_grid, ncid, "diam", unit="m")
     call bin_grid_output_netcdf(bc_grid, ncid, "bc_frac", unit="1")
     call bin_grid_output_netcdf(su_grid, ncid, "su_frac", unit="1")
     call bin_grid_output_netcdf(sc_grid, ncid, "sc", unit="1")
     call bin_grid_output_netcdf(ph_grid, ncid, "ph", unit="1")
     call bin_grid_output_netcdf(entropy_grid, ncid, "entropy", unit="1")

     call stats_1d_output_netcdf(stats_num_dist, ncid, "num_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_num_dist)

     call stats_2d_output_netcdf(stats_diam_bc_dist, ncid, "diam_bc_dist", &
          dim_name_1="diam", dim_name_2="bc_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_bc_dist)

     call stats_2d_output_netcdf(stats_diam_su_dist, ncid, "diam_su_dist", &
          dim_name_1="diam", dim_name_2="su_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_su_dist)

     call stats_2d_output_netcdf(stats_diam_su_dist_wet, ncid, "diam_su_dist_wet", &
          dim_name_1="diam", dim_name_2="su_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_su_dist_wet)

     call stats_2d_output_netcdf(stats_diam_ph_dist_acid, ncid, "diam_ph_dist_acid", &
          dim_name_1="diam", dim_name_2="ph", unit="C")
     call stats_2d_clear(stats_diam_ph_dist_acid)

     call stats_2d_output_netcdf(stats_diam_ph_dist_neut, ncid, "diam_ph_dist_neut", &
          dim_name_1="diam", dim_name_2="ph", unit="C")
     call stats_2d_clear(stats_diam_ph_dist_neut)

     call stats_2d_output_netcdf(stats_diam_sc_dist, ncid, "diam_sc_dist", &
          dim_name_1="diam", dim_name_2="sc", unit="m^{-3}")
     call stats_2d_clear(stats_diam_sc_dist)

     call stats_1d_output_netcdf(stats_entropy_dist, ncid, "entropy_dist", &
          dim_name="entropy", unit="m^{-3}")
     call stats_1d_clear(stats_entropy_dist)
     
     ! wet distribution
     call stats_1d_output_netcdf(stats_wet_num_dist, ncid, "wet_num_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_wet_num_dist)

     ! Species size distribution
     call stats_1d_output_netcdf(stats_bc_dist_dry, ncid, "bc_mass_dry_dist", &
          dim_name="diam", unit="kg /micro m")
     call stats_1d_clear(stats_bc_dist_dry)
     call stats_1d_output_netcdf(stats_bc_dist_wet, ncid, "bc_mass_wet_dist", &
          dim_name="diam", unit="kg /micro m")
     call stats_1d_clear(stats_bc_dist_wet)


     call stats_1d_output_netcdf(stats_Hp_dist_dry, ncid, "Hp_mass_dry_dist", &
          dim_name="diam", unit="kg /micro m")
     call stats_1d_clear(stats_Hp_dist_dry)
     call stats_1d_output_netcdf(stats_Hp_dist_wet, ncid, "Hp_mass_wet_dist", &
          dim_name="diam", unit="kg /micro m")
     call stats_1d_clear(stats_Hp_dist_wet)

     call stats_1d_output_netcdf(stats_h2o_dist_dry, ncid, "h2o_mass_dry_dist", &
          dim_name="diam", unit="kg /micro m")
     call stats_1d_clear(stats_h2o_dist_dry)
     call stats_1d_output_netcdf(stats_h2o_dist_wet, ncid, "h2o_mass_wet_dist", &
          dim_name="diam", unit="kg /micro m")
     call stats_1d_clear(stats_h2o_dist_wet)
   
     call stats_1d_output_netcdf(stats_su_dist_dry_acid, ncid, "su_acid_dry_dist", &
          dim_name="diam", unit="kg /micro m")
     call stats_1d_clear(stats_su_dist_dry_acid)
     
     call stats_1d_output_netcdf(stats_su_dist_dry_neut, ncid, "su_neut_dry_dist", &
          dim_name="diam", unit="kg /micro m")
     call stats_1d_clear(stats_su_dist_dry_neut)

   
     call stats_1d_output_netcdf(stats_su_dist_wet, ncid, "su_mass_wet_dist", &
          dim_name="diam", unit="kg /micro m")
     call stats_1d_clear(stats_su_dist_wet)

     call stats_1d_output_netcdf(stats_no3_dist_dry, ncid, "no3_mass_dry_dist", &
          dim_name="diam", unit="kg /micro m")
     call stats_1d_clear(stats_no3_dist_dry)
     call stats_1d_output_netcdf(stats_no3_dist_wet, ncid, "no3_mass_wet_dist", &
          dim_name="diam", unit="kg /micro m")
     call stats_1d_clear(stats_no3_dist_wet)

     call stats_1d_output_netcdf(stats_nh4_dist_dry_neut, ncid, "nh4_neut_dry_dist", &
          dim_name="diam", unit="kg /micro m")
     call stats_1d_clear(stats_nh4_dist_dry_neut)
     
     call stats_1d_output_netcdf(stats_nh4_dist_dry_acid, ncid, "nh4_acid_dry_dist",& 
          dim_name="diam", unit="kg /micro m")
     call stats_1d_clear(stats_nh4_dist_dry_acid)

     call stats_1d_output_netcdf(stats_nh4_dist_wet, ncid, "nh4_mass_wet_dist", &
          dim_name="diam", unit="kg /micro m")
     call stats_1d_clear(stats_nh4_dist_wet)
    
     call stats_1d_output_netcdf(stats_oc_dist_dry, ncid, "oc_mass_dry_dist",&
          dim_name="diam", unit="kg /micro m")
     call stats_1d_clear(stats_oc_dist_dry)
     call stats_1d_output_netcdf(stats_oc_dist_wet, ncid, "oc_mass_wet_dist",&
          dim_name="diam", unit="kg /micro m")
     call stats_1d_clear(stats_oc_dist_wet)
     
     call pmc_nc_close(ncid)
  end do

  call make_filename(out_filename, prefix, "_process.nc")
  call pmc_nc_open_write(out_filename, ncid)
  call pmc_nc_write_info(ncid, uuid, "1_urban_plume process")
  call pmc_nc_write_real_1d(ncid, times, "time", dim_name="time", unit="s")
  call stats_1d_output_netcdf(stats_tot_num_conc, ncid, "tot_num_conc", &
       dim_name="time", unit="m^{-3}")
  call stats_1d_output_netcdf(stats_tot_mass_conc, ncid,   "tot_mass_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_totdry_mass_conc, ncid,"totdry_mass_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_bc_conc, ncid, "tot_bc_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_oc_conc, ncid, "tot_oc_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_so4_conc, ncid, "tot_so4_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_no3_conc, ncid, "tot_no3_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_nh4_conc, ncid, "tot_nh4_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_h2o_conc, ncid, "tot_h2o_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_Hp_conc, ncid, "tot_Hp_conc", &
       dim_name="time", unit="kg m^{-3}")

  !%%%%%%%%%%%%SO4 RELATED SPECIES%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  call stats_1d_output_netcdf(stats_tot_aO3_conc, ncid, "tot_aO3_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_aH2O2_conc, ncid, "tot_aH2O2_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_aHO2_conc, ncid, "tot_aHO2_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_FEpp_conc, ncid, "tot_FEpp_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_FEppp_conc, ncid, "tot_FEppp_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_FEOHpp_conc, ncid, "tot_FEOHpp_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_FEC2O4p_conc, ncid, "tot_FEC2O4p_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_BRm_conc, ncid, "tot_BRm_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_aHNO4_conc, ncid, "tot_aHNO4_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_aHCHO_conc, ncid, "tot_aHCHO_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_aSO2_conc, ncid,  "tot_aSO2_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_HSO4m_conc, ncid, "tot_HSO4m_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_HSO3m_conc, ncid, "tot_HSO3m_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_SO3m_conc, ncid, "tot_SO3m_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_SO4m_conc, ncid, "tot_SO4m_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_SO3mm_conc, ncid, "tot_SO3mm_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_SO5m_conc, ncid, "tot_SO5m_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_HSO5m_conc, ncid, "tot_HSO5m_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_SO5O2Hm_conc, ncid, "tot_SO5O2Hm_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_SO5O2mm_conc, ncid, "tot_SO5O2mm_conc", &
       dim_name="time", unit="kg m^{-3}")
  !%%%%%%%%%%%%SO4 RELATED SPECIES%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !%%%%%%%%%%%%Nitric RELATED SPECIES%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  call stats_1d_output_netcdf(stats_tot_aNO3_conc, ncid, "tot_aNO3_conc",&
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_aN2O5_conc, ncid, "tot_aN2O5_conc",&
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_aNO_conc, ncid, "tot_aNO_conc",&
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_NO2p_conc, ncid, "tot_NO2p_conc",&
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_NO4m_conc, ncid, "tot_NO4m_conc",&
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_NO2m_conc, ncid, "tot_NO2m_conc",&
       dim_name="time", unit="kg m^{-3}")
 !%%%%%%%%%%%%Nitric RELATED SPECIES%%%%%%%%
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 !%%%%%%%%%%%%TEST%%%%%%%%%%%%%%%%%%%%%%%
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 call stats_1d_output_netcdf(stats_tot_aH2C2O4_conc, ncid, "tot_aH2C2O4_conc",&
       dim_name="time", unit="kg m^{-3}")
 call stats_1d_output_netcdf(stats_tot_C2O4m_conc, ncid,   "tot_C2O4m_conc",&
       dim_name="time", unit="kg m^{-3}")
 call stats_1d_output_netcdf(stats_tot_HC2O4m_conc, ncid, "tot_HC2O4m_conc",&
       dim_name="time", unit="kg m^{-3}")
 call stats_1d_output_netcdf(stats_tot_aCO2_conc, ncid, "tot_aCO2_conc",&
       dim_name="time", unit="kg m^{-3}")
 call stats_1d_output_netcdf(stats_tot_aO2_conc, ncid, "tot_aO2_conc",&
       dim_name="time", unit="kg m^{-3}")
 !%%%%%%%%%%%%Nitric RELATED SPECIES%%%%%
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  call stats_1d_output_netcdf(stats_tot_aNO2_conc, ncid, "tot_aNO2_conc",&
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_aOH_conc, ncid, "tot_aOH_conc",&
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_aNH3_conc, ncid, "tot_aNH3_conc",&
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_entropy, ncid, "tot_entropy", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_rh,  ncid, "rh", &
       dim_name="time", unit="%")
  call stats_1d_output_netcdf(stats_pre, ncid, "pre", &
       dim_name="time", unit="pa")
  call stats_1d_output_netcdf(stats_temp, ncid, "temp", &
       dim_name="time", unit="K")
  call stats_1d_output_netcdf(stats_den, ncid, "density", &
       dim_name="time", unit="kg m^-3")
  call stats_1d_output_netcdf(stats_spres, ncid, "satu_pres", &
       dim_name="time", unit="pa")
  call stats_1d_output_netcdf(stats_tot_entropy_averaged, ncid, &
       "tot_entropy_averaged", dim_name="time", unit="m^{-3}")
  call pmc_nc_close(ncid)

  call bin_grid_allocate(diam_grid)
  call aero_data_deallocate(aero_data)
  call aero_state_deallocate(aero_state)
  call env_state_deallocate(env_state)

  call pmc_mpi_finalize()

end program process
