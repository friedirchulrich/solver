module variablen


!! Variablen zur Image-Kontrolle
INTEGER::M,N,O,m_kop,n_kop,o_kop,myrank,numprocs,m_pos,n_pos,o_pos,R,S,T_p&
,m_local,n_local,o_local,image_m_1,image_n_1,image_o_1,image_m_max,image_n_max,&
image_o_max
INTEGER::sync_m_up,sync_m_down,sync_n_up,sync_n_down,sync_o_up,sync_o_down,nlocal
INTEGER,allocatable::sync_m_array(:),sync_n_array(:),sync_o_array(:)

!! Ausgabe-parameter
INTEGER,PARAMETER::prtrows=10
INTEGER,PARAMETER::prtcols=10
LOGICAL::start

!! Felder der Physikgrößen
!! Co-arrays
REAL,allocatable::rho(:,:,:)[:,:,:],rho_u(:,:,:)[:,:,:],rho_E(:,:,:)[:,:,:],rho_v(:,:,:)[:,:,:]&
,E(:,:,:)[:,:,:],p(:,:,:)[:,:,:],Epu(:,:,:)[:,:,:],Epv(:,:,:)[:,:,:],Epw(:,:,:)[:,:,:],&
E_k(:,:,:)[:,:,:],vis(:,:,:)[:,:,:],mu(:,:,:)[:,:,:]&
,a(:,:,:)[:,:,:],rho_w(:,:,:)[:,:,:],k(:,:,:)[:,:,:],&
rho_neu(:,:,:)[:,:,:],rho_u_neu(:,:,:)[:,:,:],rho_E_neu(:,:,:)[:,:,:],rho_v_neu(:,:,:)[:,:,:],&
rho_w_neu(:,:,:)[:,:,:]
REAL,allocatable::T(:,:,:)[:,:,:],u(:,:,:)[:,:,:],v(:,:,:)[:,:,:],DDX_DDY(:,:,:)[:,:,:],w(:,:,:)[:,:,:]
REAL,allocatable::dphi(:,:,:)[:,:,:],solve(:,:,:)[:,:,:],&
dphi1(:,:,:)[:,:,:],dphi2(:,:,:)[:,:,:],dphi3(:,:,:)[:,:,:],dphi4(:,:,:)[:,:,:]
REAL,allocatable::du(:,:,:)[:,:,:],dv(:,:,:)[:,:,:]
REAL,allocatable::poisson(:,:,:)[:,:,:]
REAL,allocatable::uinit(:,:,:)[:,:,:],rinit(:,:,:)[:,:,:],vinit(:,:,:)[:,:,:],tinit(:,:,:)[:,:,:]
REAL,allocatable::drudx(:,:,:)[:,:,:],drvdy(:,:,:)[:,:,:],&
dudy(:,:,:)[:,:,:],dvdy(:,:,:)[:,:,:],dudx(:,:,:)[:,:,:],dvdx(:,:,:)[:,:,:],&
rho_bc(:,:,:)[:,:,:],T_bc(:,:,:)[:,:,:],u_bc(:,:,:)[:,:,:],v_bc(:,:,:)[:,:,:]

!! Felder für Gleichungslösung
!in x
REAL,allocatable::b_x(:),c_x(:),a_x(:)[:,:,:],f_x(:)[:,:,:],g_x(:)[:,:,:],x_x(:,:,:)[:,:,:],y_x(:,:,:)[:,:,:]
REAL,allocatable::a_x_s(:)[:,:,:],f_x_s(:)[:,:,:],g_x_s(:)[:,:,:],y_x_s(:,:,:)[:,:,:],x_x_s(:,:,:)[:,:,:]
REAL,allocatable::alpha_x(:,:)[:,:,:],beta_x(:)[:,:,:],z_x(:)
!in y
REAL,allocatable::b_y(:),c_y(:),a_y(:)[:,:,:],f_y(:)[:,:,:],g_y(:)[:,:,:],x_y(:,:,:)[:,:,:],y_y(:,:,:)[:,:,:]
REAL,allocatable::a_y_s(:)[:,:,:],f_y_s(:)[:,:,:],g_y_s(:)[:,:,:],y_y_s(:,:,:)[:,:,:],x_y_s(:,:,:)[:,:,:]
REAL,allocatable::alpha_y(:,:)[:,:,:],beta_y(:)[:,:,:],z_y(:)
!in z
REAL,allocatable::b_z(:),c_z(:),a_z(:)[:,:,:],f_z(:)[:,:,:],g_z(:)[:,:,:],x_z(:,:,:)[:,:,:],y_z(:,:,:)[:,:,:]
REAL,allocatable::a_z_s(:)[:,:,:],f_z_s(:)[:,:,:],g_z_s(:)[:,:,:],y_z_s(:,:,:)[:,:,:],x_z_s(:,:,:)[:,:,:]
REAL,allocatable::alpha_z(:,:)[:,:,:],beta_z(:)[:,:,:],z_z(:)

!! Koeffizienten für kompakte Differenzen
REAL::alpha_i,beta_i,gamma_i,alpha_nr,beta_nr,gamma_nr,beta_r,gamma_r
REAL::q_i,a_i,b_i,c_i,d_i,e_i,q_nr,b_nr,c_nr,d_nr,q_r,c_r,d_r,e_r
REAL::alpha2_i,beta2_i,gamma2_i,alpha2_nr,beta2_nr,gamma2_nr,beta2_r,gamma2_r
REAL::q2_i,a2_i,b2_i,c2_i,d2_i,e2_i,q2_nr,b2_nr,c2_nr,d2_nr,q2_r,c2_r,d2_r,e2_r
REAL::qxi_i,qxi_nr,qxi_r,qyi_i,qyi_nr,qyi_r,qzi_i,qzi_nr,qzi_r
REAL::qxi2_i,qxi2_nr,qxi2_r,qyi2_i,qyi2_nr,qyi2_r,qzi2_i,qzi2_nr,qzi2_r

!! Physik
REAL::dx,dy,dz,dt,Ma
REAL::R_gas,cv_gas,kappa,cv,cp, T_inf,Pr,gamma
REAL::u_0,rho_0,p_0,E_0,vis0
REAL::mu0,T0,Suth !mu0, T0 , Sutherland constant
!!
REAL::total_energy[1:*],total_konti[1:*],E_k_total[1:*]
REAL,allocatable::kontinuitaet(:,:)[:]

!!Hilfsvariablen
REAL,allocatable::E_plus_p(:,:,:)[:,:,:]

!!spannungstensor
REAL,allocatable:: tau_xx(:,:,:)[:,:,:],tau_yy(:,:,:)[:,:,:],tau_zz(:,:,:)[:,:,:]
REAL,allocatable:: tau_xy(:,:,:)[:,:,:],tau_xz(:,:,:)[:,:,:],tau_zy(:,:,:)[:,:,:]
REAL,allocatable:: delta_u(:,:,:)[:,:,:],delta_v(:,:,:)[:,:,:],delta_w(:,:,:)[:,:,:]
end module
