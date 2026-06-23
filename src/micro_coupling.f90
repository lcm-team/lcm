SUBROUTINE micro_coupling
    !
    !-- Prepare LES fields (u, v, w, tabs) to be used in the LCM
        USE micro_vars,         ONLY: u_LCM, v_LCM, w_LCM, tabs_LCM, rho_air,    &
                                      pres_LCM, dtn_LCM, nx_LCM, ny_LCM, nz_LCM, &
                                      x1u, x2u, y1u, y2u, ux1, ux2, uy1, uy2,    &
                                      x1v, x2v, y1v, y2v, vx1, vx2, vy1, vy2,    &
                                      x1w, x2w, y1w, y2w, wx1, wx2, wy1, wy2
    
        USE grid,               ONLY: dtn, nx, ny, nzm
        USE vars,               ONLY: u, v, w, t, rho, pres, tabs, gamaz
    
        USE microphysics,       ONLY: micro_field, iTlcm, iTllcm
    
        IMPLICIT NONE
    
        INTEGER :: k
    
        CALL t_startf ('micro_coupling')
    
        u_LCM(1:nx_LCM,1:ny_LCM,1:nz_LCM) = u(1:nx,1:ny,1:nzm)
        v_LCM(1:nx_LCM,1:ny_LCM,1:nz_LCM) = v(1:nx,1:ny,1:nzm)
        w_LCM(1:nx_LCM,1:ny_LCM,1:nz_LCM) = w(1:nx,1:ny,1:nzm)
    
        CALL task_exchange_LCM( u_LCM, x1u, x2u, y1u, y2u, nz_LCM, ux1, ux2, uy1, uy2, 1 )
        CALL task_exchange_LCM( v_LCM, x1v, x2v, y1v, y2v, nz_LCM, vx1, vx2, vy1, vy2, 2 )
        CALL task_exchange_LCM( w_LCM, x1w, x2w, y1w, y2w, nz_LCM, wx1, wx2, wy1, wy2, 3 )
    !
    !-- Get effect of surface fluxes, large-scale forcings, subsidence, radiation, etc.
    !-- from SAM's prognostic equation for t (liquid water static energy)
        micro_field(:,:,:,iTlcm) = micro_field(:,:,:,iTlcm) + ( t(:,:,:) - micro_field(:,:,:,iTllcm) )
        DO k = 1, nz_LCM
           tabs_LCM(1:nx_LCM,1:ny_LCM,k) = micro_field(1:nx_LCM,1:ny_LCM,k,iTlcm) - gamaz(k)
        END DO
    
        rho_air(1:nz_LCM)  = rho(1:nzm)
        pres_LCM(1:nz_LCM) = pres(1:nzm)
    
        dtn_LCM = dtn
    
        CALL t_stopf ('micro_coupling')
    
     END SUBROUTINE micro_coupling
    